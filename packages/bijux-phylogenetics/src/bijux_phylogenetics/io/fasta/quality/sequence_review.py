from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    DuplicateSequencePolicyAction,
    DuplicateSequencePolicyReport,
    SequenceQualityRankingReport,
    SequenceQualityRankingRow,
)

from ..cleaning import (
    detect_identical_duplicate_sequences,
    detect_near_duplicate_sequences,
)
from ..records import summarise_fasta


def build_duplicate_sequence_policy_report(
    path: Path,
    *,
    near_duplicate_threshold: float = 0.99,
) -> DuplicateSequencePolicyReport:
    """Build reviewer-facing policy recommendations for duplicate sequences."""
    exact_duplicates = detect_identical_duplicate_sequences(path)
    near_duplicates = detect_near_duplicate_sequences(
        path, identity_threshold=near_duplicate_threshold
    )
    actions: list[DuplicateSequencePolicyAction] = []
    warnings: list[str] = []
    if exact_duplicates:
        warnings.append(
            "exact duplicate sequences should be deduplicated or explicitly justified before inference"
        )
        for group in exact_duplicates:
            actions.append(
                DuplicateSequencePolicyAction(
                    action="collapse_exact_duplicates",
                    rationale=(
                        f"retain one representative such as {group.identifiers[0]} unless metadata show that "
                        "the duplicated labels represent distinct biological samples"
                    ),
                    affected_identifiers=group.identifiers,
                )
            )
    if near_duplicates:
        warnings.append(
            "near-duplicate sequences should be checked for replicate samples, contamination, or oversampling bias"
        )
        seen_pairs: set[tuple[str, ...]] = set()
        for pair in near_duplicates:
            identifiers = tuple(sorted((pair.left_identifier, pair.right_identifier)))
            if identifiers in seen_pairs:
                continue
            seen_pairs.add(identifiers)
            actions.append(
                DuplicateSequencePolicyAction(
                    action="review_near_duplicates",
                    rationale="inspect metadata and sampling provenance before keeping highly similar sequences together in inference",
                    affected_identifiers=list(identifiers),
                )
            )
    if not warnings:
        warnings.append("no duplicate-sequence policy actions are currently required")
    return DuplicateSequencePolicyReport(
        path=path,
        exact_duplicate_groups=exact_duplicates,
        near_duplicate_pairs=near_duplicates,
        policy_actions=actions,
        warnings=warnings,
    )


def build_sequence_quality_ranking(path: Path) -> SequenceQualityRankingReport:
    """Rank aligned sequences by transparent quality burdens."""
    summary = summarise_fasta(path)
    composition_outlier_ids = {row.identifier for row in summary.composition_outliers}
    exact_duplicate_ids = {
        identifier
        for group in summary.duplicate_sequence_groups
        for identifier in group.identifiers
    }
    near_duplicate_ids = {
        identifier
        for pair in summary.near_duplicate_pairs
        for identifier in (pair.left_identifier, pair.right_identifier)
    }
    uncertainty_by_id = {
        row.identifier: row for row in summary.per_sequence_uncertainty
    }
    ranked: list[tuple[float, str, SequenceQualityRankingRow]] = []
    for identifier in summary.ids:
        uncertainty = uncertainty_by_id[identifier]
        composition_outlier = identifier in composition_outlier_ids
        if identifier in exact_duplicate_ids:
            duplicate_status = "exact_duplicate"
        elif identifier in near_duplicate_ids:
            duplicate_status = "near_duplicate"
        else:
            duplicate_status = "unique"
        penalty = (
            uncertainty.missing_fraction * 40.0
            + uncertainty.gap_fraction * 25.0
            + uncertainty.ambiguity_fraction * 20.0
            + (10.0 if composition_outlier else 0.0)
            + (
                10.0
                if duplicate_status == "exact_duplicate"
                else 5.0
                if duplicate_status == "near_duplicate"
                else 0.0
            )
        )
        score = round(max(0.0, 100.0 - penalty), 3)
        note_parts: list[str] = []
        if uncertainty.missing_fraction > 0.0:
            note_parts.append("missing data")
        if uncertainty.gap_fraction > 0.0:
            note_parts.append("gaps")
        if uncertainty.ambiguity_fraction > 0.0:
            note_parts.append("ambiguity codes")
        if composition_outlier:
            note_parts.append("composition outlier")
        if duplicate_status != "unique":
            note_parts.append(duplicate_status.replace("_", " "))
        note = (
            "quality burdens: " + ", ".join(note_parts)
            if note_parts
            else "no major quality burdens detected"
        )
        ranked.append(
            (
                score,
                identifier,
                SequenceQualityRankingRow(
                    identifier=identifier,
                    rank=0,
                    score=score,
                    missing_fraction=uncertainty.missing_fraction,
                    gap_fraction=uncertainty.gap_fraction,
                    ambiguity_fraction=uncertainty.ambiguity_fraction,
                    composition_outlier=composition_outlier,
                    duplicate_status=duplicate_status,
                    note=note,
                ),
            )
        )
    ranked.sort(key=lambda item: (item[0], item[1]))
    rows = [
        SequenceQualityRankingRow(
            identifier=row.identifier,
            rank=index,
            score=row.score,
            missing_fraction=row.missing_fraction,
            gap_fraction=row.gap_fraction,
            ambiguity_fraction=row.ambiguity_fraction,
            composition_outlier=row.composition_outlier,
            duplicate_status=row.duplicate_status,
            note=row.note,
        )
        for index, (_, _, row) in enumerate(ranked, start=1)
    ]
    warnings = (
        ["lower-ranked sequences should be reviewed before publication or inference"]
        if rows and any(row.score < 85.0 for row in rows)
        else []
    )
    if not summary.near_duplicate_scan_performed:
        warnings.append(
            "near-duplicate sequence ranking was skipped because the alignment exceeds the governed pairwise review threshold"
        )
    return SequenceQualityRankingReport(path=path, rows=rows, warnings=warnings)
