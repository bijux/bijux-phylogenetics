from __future__ import annotations

import csv
from pathlib import Path

from ..tree_sets.clade_support import (
    _support_classification,
    compute_clade_frequency_table,
)
from ..tree_sets.inventory import (
    _analyze_tree_set,
    _require_exact_taxa,
    _require_tree_set,
    _TreeSetAnalysis,
    _validate_same_taxa,
)
from ..tree_sets.topology import (
    _clade_counts,
    _clade_signature,
    _clades_conflict,
    _format_clade,
)
from .models import (
    CladeCredibilityConflict,
    CladeCredibilityConflictReport,
    TaxonPlacementSignature,
    UncertaintyAwareCladeConclusion,
    UncertaintyAwareConclusionSummaryReport,
    UnstableClade,
    UnstableCladeReport,
    UnstableTaxaReport,
    UnstableTaxon,
)


def _build_unstable_clade_report(analysis: _TreeSetAnalysis) -> UnstableCladeReport:
    _require_exact_taxa(analysis)
    counts = analysis.clade_counts or {}
    all_clades = set(counts)
    tree_count = len(analysis.trees)
    unstable_clades = [
        UnstableClade(
            clade=_format_clade(clade),
            tree_count=count,
            frequency=round(count / tree_count, 15),
            conflict_count=len(
                conflicts := sorted(
                    _format_clade(other)
                    for other in all_clades
                    if _clades_conflict(clade, other)
                )
            ),
            instability_score=round(
                min(count / tree_count, 1.0 - (count / tree_count)), 15
            ),
            support_classification=_support_classification(
                round(count / tree_count, 15),
                len(conflicts),
            ),
            conflicting_clades=conflicts,
        )
        for clade, count in sorted(
            counts.items(), key=lambda item: _format_clade(item[0])
        )
        if count < tree_count
    ]
    unstable_clades.sort(
        key=lambda row: (-row.instability_score, -row.conflict_count, row.clade)
    )
    return UnstableCladeReport(
        path=analysis.path,
        tree_count=tree_count,
        processing=analysis.processing,
        clades=unstable_clades,
    )


def detect_unstable_taxa(path: Path) -> UnstableTaxaReport:
    """Report taxa whose placement signatures vary across trees in a set."""
    analysis = _analyze_tree_set(path)
    trees = analysis.trees
    shared_taxa = set(_require_exact_taxa(analysis))
    taxa: list[UnstableTaxon] = []
    for taxon in sorted(shared_taxa):
        signature_counts: dict[str, int] = {}
        for tree in trees:
            signature = _clade_signature(tree, shared_taxa, taxon)
            signature_counts[signature] = signature_counts.get(signature, 0) + 1
        if len(signature_counts) < 2:
            continue
        placements = [
            TaxonPlacementSignature(
                signature=signature,
                tree_count=count,
                frequency=round(count / len(trees), 15),
            )
            for signature, count in sorted(
                signature_counts.items(),
                key=lambda item: (-item[1], item[0]),
            )
        ]
        taxa.append(
            UnstableTaxon(
                taxon=taxon,
                unique_placements=len(signature_counts),
                dominant_frequency=placements[0].frequency,
                instability_score=round(1.0 - placements[0].frequency, 15),
                placements=placements,
            )
        )
    taxa.sort(
        key=lambda row: (-row.instability_score, -row.unique_placements, row.taxon)
    )
    return UnstableTaxaReport(
        path=path,
        tree_count=len(trees),
        processing=analysis.processing,
        taxa=taxa,
    )


def detect_unstable_clades(path: Path) -> UnstableCladeReport:
    """Report non-unanimous clades and their conflicting alternatives."""
    return _build_unstable_clade_report(_analyze_tree_set(path))


def write_unstable_clade_table(path: Path, report: UnstableCladeReport) -> Path:
    """Write unstable clades and their conflicting alternatives as a TSV table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "clade",
                "tree_count",
                "frequency",
                "conflict_count",
                "instability_score",
                "support_classification",
                "conflicting_clades",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.clades:
            writer.writerow(
                {
                    "clade": row.clade,
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                    "conflict_count": row.conflict_count,
                    "instability_score": format(row.instability_score, ".15g"),
                    "support_classification": row.support_classification,
                    "conflicting_clades": ",".join(row.conflicting_clades),
                }
            )
    return path


def summarize_clade_credibility_conflicts(
    path: Path,
    *,
    credibility_threshold: float = 0.5,
) -> CladeCredibilityConflictReport:
    """Identify mutually incompatible clades that both achieve high posterior credibility."""
    if not 0.0 < credibility_threshold < 1.0:
        raise ValueError(
            f"credibility_threshold must be between 0 and 1, got {credibility_threshold}"
        )
    _, trees = _require_tree_set(path)
    shared_taxa = set(_validate_same_taxa(trees))
    counts = _clade_counts(trees, shared_taxa)
    frequencies = {
        clade: round(count / len(trees), 15) for clade, count in counts.items()
    }
    high_credibility = sorted(
        [
            clade
            for clade, frequency in frequencies.items()
            if frequency >= credibility_threshold
        ],
        key=_format_clade,
    )
    conflicts: list[CladeCredibilityConflict] = []
    for index, left_clade in enumerate(high_credibility):
        for right_clade in high_credibility[index + 1 :]:
            if not _clades_conflict(left_clade, right_clade):
                continue
            conflicts.append(
                CladeCredibilityConflict(
                    left_clade=_format_clade(left_clade),
                    left_frequency=frequencies[left_clade],
                    right_clade=_format_clade(right_clade),
                    right_frequency=frequencies[right_clade],
                    combined_frequency=round(
                        frequencies[left_clade] + frequencies[right_clade], 15
                    ),
                )
            )
    conflicts.sort(
        key=lambda row: (-row.combined_frequency, row.left_clade, row.right_clade)
    )
    return CladeCredibilityConflictReport(
        path=path,
        tree_count=len(trees),
        credibility_threshold=credibility_threshold,
        high_credibility_clade_count=len(high_credibility),
        conflict_count=len(conflicts),
        conflicts=conflicts,
    )


def write_clade_credibility_conflict_table(
    path: Path,
    report: CladeCredibilityConflictReport,
) -> Path:
    """Write high-credibility clade conflicts as a TSV table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "left_clade",
                "left_frequency",
                "right_clade",
                "right_frequency",
                "combined_frequency",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.conflicts:
            writer.writerow(
                {
                    "left_clade": row.left_clade,
                    "left_frequency": format(row.left_frequency, ".15g"),
                    "right_clade": row.right_clade,
                    "right_frequency": format(row.right_frequency, ".15g"),
                    "combined_frequency": format(row.combined_frequency, ".15g"),
                }
            )
    return path


def summarize_uncertainty_aware_conclusions(
    path: Path,
    *,
    robust_threshold: float = 0.9,
    uncertain_min_frequency: float = 0.3,
    uncertain_max_frequency: float = 0.7,
    credibility_threshold: float = 0.5,
) -> UncertaintyAwareConclusionSummaryReport:
    """Classify clade-level conclusions as robust, uncertain, or conflict-prone."""
    if not 0.0 < robust_threshold <= 1.0:
        raise ValueError(
            f"robust_threshold must be between 0 and 1, got {robust_threshold}"
        )
    detect_unstable_clades(path)
    conflict_report = summarize_clade_credibility_conflicts(
        path, credibility_threshold=credibility_threshold
    )
    # Flag one ordered clade representative per incompatibility edge so the
    # conclusion summary does not double-count both sides of every conflict row.
    conflict_clades = {row.left_clade for row in conflict_report.conflicts}
    frequency_report = compute_clade_frequency_table(path)
    robust_clades: list[UncertaintyAwareCladeConclusion] = []
    uncertain_clades: list[UncertaintyAwareCladeConclusion] = []
    conflicting_clades: list[UncertaintyAwareCladeConclusion] = []
    for row in frequency_report.clade_frequencies:
        if row.clade in conflict_clades:
            conflicting_clades.append(
                UncertaintyAwareCladeConclusion(
                    clade=row.clade,
                    frequency=row.frequency,
                    conclusion="conflict-prone",
                    rationale="clade reaches high posterior frequency but is incompatible with another high-credibility clade",
                )
            )
            continue
        if row.frequency >= robust_threshold:
            robust_clades.append(
                UncertaintyAwareCladeConclusion(
                    clade=row.clade,
                    frequency=row.frequency,
                    conclusion="robust",
                    rationale="clade remains near-fixed across the posterior tree set",
                )
            )
            continue
        if uncertain_min_frequency <= row.frequency <= uncertain_max_frequency:
            uncertain_clades.append(
                UncertaintyAwareCladeConclusion(
                    clade=row.clade,
                    frequency=row.frequency,
                    conclusion="uncertain",
                    rationale="clade holds intermediate support and should not anchor strong biological interpretation",
                )
            )
    robust_clades.sort(key=lambda row: (-row.frequency, row.clade))
    uncertain_clades.sort(key=lambda row: (-row.frequency, row.clade))
    conflicting_clades.sort(key=lambda row: (-row.frequency, row.clade))
    return UncertaintyAwareConclusionSummaryReport(
        path=path,
        tree_count=frequency_report.tree_count,
        robust_clade_count=len(robust_clades),
        uncertain_clade_count=len(uncertain_clades),
        conflicting_clade_count=len(conflicting_clades),
        robust_clades=robust_clades,
        uncertain_clades=uncertain_clades,
        conflicting_clades=conflicting_clades,
    )


def write_uncertainty_conclusion_table(
    path: Path,
    report: UncertaintyAwareConclusionSummaryReport,
) -> Path:
    """Write robust, uncertain, and conflict-prone clade conclusions as a TSV table."""
    rows = [
        *report.robust_clades,
        *report.uncertain_clades,
        *report.conflicting_clades,
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["clade", "frequency", "conclusion", "rationale"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "clade": row.clade,
                    "frequency": format(row.frequency, ".15g"),
                    "conclusion": row.conclusion,
                    "rationale": row.rationale,
                }
            )
    return path
