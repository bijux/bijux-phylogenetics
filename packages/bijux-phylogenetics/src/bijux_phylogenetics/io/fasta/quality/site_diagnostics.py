from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    AlignmentAmbiguousColumnReport,
    AlignmentLowInformationReport,
    AlignmentMissingDataConcentration,
    AlignmentSummary,
    AmbiguousAlignmentColumn,
)

from ..core import _validate_fraction_threshold
from ..records import summarise_fasta

LOW_INFORMATION_SITE_THRESHOLD = 1
LOW_INFORMATION_FRACTION_THRESHOLD = 0.0


def assess_alignment_low_information_from_summary(
    summary: AlignmentSummary,
    *,
    minimum_informative_sites: int,
    minimum_informative_fraction: float,
) -> AlignmentLowInformationReport:
    informative_fraction = (
        0.0
        if summary.alignment_length == 0
        else round(
            summary.parsimony_informative_site_count / summary.alignment_length,
            15,
        )
    )
    reasons: list[str] = []
    if (
        summary.parsimony_informative_site_count < minimum_informative_sites
        and summary.variable_site_count == 0
    ):
        reasons.append(
            "alignment has fewer parsimony-informative sites than the minimum threshold for defensible inference"
        )
    if (
        summary.variable_site_count == 0
        and informative_fraction < minimum_informative_fraction
    ):
        reasons.append(
            "alignment has a very low parsimony-informative-site fraction and may not support stable inference"
        )
    return AlignmentLowInformationReport(
        sequence_count=summary.sequence_count,
        alignment_length=summary.alignment_length,
        parsimony_informative_site_count=summary.parsimony_informative_site_count,
        parsimony_informative_fraction=informative_fraction,
        threshold_site_count=minimum_informative_sites,
        threshold_fraction=minimum_informative_fraction,
        low_information=bool(reasons),
        reasons=reasons,
    )


def build_ambiguous_alignment_column_report_from_summary(
    path: Path,
    summary: AlignmentSummary,
    *,
    threshold: float,
) -> AlignmentAmbiguousColumnReport:
    uncertainty_by_position = {
        row.position: row for row in summary.per_site_uncertainty
    }
    rows: list[AmbiguousAlignmentColumn] = []
    for position in range(1, summary.alignment_length + 1):
        uncertainty = uncertainty_by_position[position]
        ambiguity_burden = (
            uncertainty.ambiguity_fraction
            + uncertainty.missing_fraction
            + uncertainty.gap_fraction
        )
        if ambiguity_burden < threshold:
            continue
        rows.append(
            AmbiguousAlignmentColumn(
                position=position,
                ambiguity_fraction=uncertainty.ambiguity_fraction,
                missing_fraction=uncertainty.missing_fraction,
                gap_fraction=uncertainty.gap_fraction,
                comparable_fraction=round(max(0.0, 1.0 - ambiguity_burden), 15),
                note="column is dominated by ambiguity, missing data, or gaps and should be reviewed for masking",
            )
        )
    warnings = (
        [
            "alignment contains ambiguity-heavy columns that may be unsuitable for inference without masking"
        ]
        if rows
        else []
    )
    return AlignmentAmbiguousColumnReport(
        path=path,
        threshold=threshold,
        rows=rows,
        warnings=warnings,
    )


def alignment_quality_components(summary: AlignmentSummary) -> dict[str, float]:
    informative_density = (
        0.0
        if summary.alignment_length == 0
        else summary.parsimony_informative_site_count / summary.alignment_length
    )
    return {
        "missingness": round(
            max(0.0, 1.0 - min(summary.missing_data_fraction, 1.0)),
            15,
        ),
        "gap_burden": round(max(0.0, 1.0 - min(summary.gap_fraction, 1.0)), 15),
        "composition_outliers": round(
            max(
                0.0,
                1.0
                - min(
                    len(summary.composition_outliers) / max(summary.sequence_count, 1),
                    1.0,
                ),
            ),
            15,
        ),
        "duplicates": round(
            max(
                0.0,
                1.0
                - min(
                    len(summary.duplicate_sequence_groups)
                    / max(summary.sequence_count, 1),
                    1.0,
                ),
            ),
            15,
        ),
        "informative_density": round(min(informative_density / 0.2, 1.0), 15),
    }


def alignment_quality_score(components: dict[str, float]) -> float:
    return round(sum(components.values()) / max(len(components), 1) * 100.0, 3)


def summarize_missing_data_concentration(
    summary: AlignmentSummary,
    *,
    threshold: float = 0.5,
) -> AlignmentMissingDataConcentration:
    concentrated_columns = [
        row.position
        for row in summary.per_site_missingness
        if row.missing_fraction >= threshold
    ]
    longest_run = 0
    longest_start: int | None = None
    longest_end: int | None = None
    current_run = 0
    current_start: int | None = None
    previous_position: int | None = None
    for position in concentrated_columns:
        if previous_position is None or position != previous_position + 1:
            current_run = 1
            current_start = position
        else:
            current_run += 1
        if current_run > longest_run:
            longest_run = current_run
            longest_start = current_start
            longest_end = position
        previous_position = position
    maximum_missing_fraction = max(
        (row.missing_fraction for row in summary.per_site_missingness),
        default=0.0,
    )
    maximum_missing_positions = [
        row.position
        for row in summary.per_site_missingness
        if row.missing_fraction == maximum_missing_fraction
    ]
    return AlignmentMissingDataConcentration(
        threshold=threshold,
        concentrated_column_count=len(concentrated_columns),
        concentrated_column_fraction=(
            0.0
            if summary.alignment_length == 0
            else len(concentrated_columns) / summary.alignment_length
        ),
        longest_concentrated_run=longest_run,
        longest_concentrated_run_start=longest_start,
        longest_concentrated_run_end=longest_end,
        maximum_missing_fraction=maximum_missing_fraction,
        maximum_missing_positions=maximum_missing_positions,
    )


def alignment_suspicion_reasons(
    *,
    low_information: AlignmentLowInformationReport,
    missing_data_concentration: AlignmentMissingDataConcentration,
    ambiguous_column_count: int,
    over_aligned_count: int,
    under_aligned_count: int,
    invalid_character_count: int,
) -> list[str]:
    reasons: list[str] = []
    if low_information.low_information:
        reasons.append("alignment has low information content for defensible inference")
    if missing_data_concentration.longest_concentrated_run >= 2:
        reasons.append("alignment concentrates missing data into adjacent columns")
    elif missing_data_concentration.concentrated_column_count > 0:
        reasons.append("alignment contains one or more highly missing columns")
    if ambiguous_column_count > 0:
        reasons.append("alignment contains ambiguity-heavy columns")
    if over_aligned_count > 0:
        reasons.append("alignment contains suspiciously over-aligned windows")
    if under_aligned_count > 0:
        reasons.append("alignment contains suspiciously under-aligned windows")
    if invalid_character_count > 0:
        reasons.append(
            "alignment contains invalid characters for the inferred alphabet"
        )
    return reasons


def assess_alignment_low_information(
    path: Path,
    *,
    minimum_informative_sites: int = LOW_INFORMATION_SITE_THRESHOLD,
    minimum_informative_fraction: float = LOW_INFORMATION_FRACTION_THRESHOLD,
) -> AlignmentLowInformationReport:
    """Assess whether an alignment carries enough signal for defensible inference."""
    _validate_fraction_threshold(minimum_informative_fraction)
    if minimum_informative_sites < 0:
        raise ValueError(
            f"minimum_informative_sites must be non-negative, got {minimum_informative_sites}"
        )
    summary = summarise_fasta(path)
    return assess_alignment_low_information_from_summary(
        summary,
        minimum_informative_sites=minimum_informative_sites,
        minimum_informative_fraction=minimum_informative_fraction,
    )


def build_ambiguous_alignment_column_report(
    path: Path,
    *,
    threshold: float = 0.5,
) -> AlignmentAmbiguousColumnReport:
    """Report columns dominated by ambiguity, explicit missingness, or gaps."""
    _validate_fraction_threshold(threshold)
    summary = summarise_fasta(path)
    return build_ambiguous_alignment_column_report_from_summary(
        path,
        summary,
        threshold=threshold,
    )
