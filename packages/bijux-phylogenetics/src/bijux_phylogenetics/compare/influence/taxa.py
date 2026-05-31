from __future__ import annotations

from copy import deepcopy
import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation.structure import _load_tree

from ..topology.comparison import _build_tree_comparison_report
from ..topology.models import (
    SupportComparisonReport,
    TreeComparisonReport,
)
from ..topology.support import (
    _build_support_comparison_report,
)


@dataclass(slots=True)
class TaxonInfluenceRow:
    taxon: str
    retained_taxa: list[str]
    baseline_topology_equal: bool
    leave_one_out_topology_equal: bool
    baseline_same_unrooted_topology: bool
    leave_one_out_same_unrooted_topology: bool
    baseline_same_taxa_different_rooting: bool
    leave_one_out_same_taxa_different_rooting: bool
    baseline_rooted_robinson_foulds_distance: int
    leave_one_out_rooted_robinson_foulds_distance: int
    rooted_robinson_foulds_delta: int
    baseline_rooted_normalized_robinson_foulds: float
    leave_one_out_rooted_normalized_robinson_foulds: float
    rooted_normalized_robinson_foulds_delta: float
    baseline_unrooted_robinson_foulds_distance: int
    leave_one_out_unrooted_robinson_foulds_distance: int
    unrooted_robinson_foulds_delta: int
    baseline_unrooted_normalized_robinson_foulds: float
    leave_one_out_unrooted_normalized_robinson_foulds: float
    unrooted_normalized_robinson_foulds_delta: float
    baseline_support_disagreements: int
    leave_one_out_support_disagreements: int
    support_disagreement_delta: int
    baseline_conflicting_clades: int
    leave_one_out_conflicting_clades: int
    conflicting_clade_delta: int
    baseline_high_support_conflicts: int
    leave_one_out_high_support_conflicts: int
    high_support_conflict_delta: int
    topology_changed: bool
    support_changed: bool
    influence_score: float
    influence_rank: int


@dataclass(slots=True)
class TaxonInfluenceReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    baseline_topology: TreeComparisonReport
    baseline_support: SupportComparisonReport
    rows: list[TaxonInfluenceRow]


def _count_support_disagreements(report: SupportComparisonReport) -> int:
    return sum(1 for row in report.shared_clades if row.support_disagreement)


def _count_high_support_conflicts(report: SupportComparisonReport) -> int:
    return sum(
        1
        for row in report.conflicting_clades
        if row.conflict_classification == "high_support_conflict"
    )


def _build_taxon_influence_row(
    taxon: str,
    retained_taxa: list[str],
    baseline_topology: TreeComparisonReport,
    leave_one_out_topology: TreeComparisonReport,
    baseline_support: SupportComparisonReport,
    leave_one_out_support: SupportComparisonReport,
) -> TaxonInfluenceRow:
    baseline_support_disagreements = _count_support_disagreements(baseline_support)
    leave_one_out_support_disagreements = _count_support_disagreements(
        leave_one_out_support
    )
    baseline_conflicting_clades = len(baseline_support.conflicting_clades)
    leave_one_out_conflicting_clades = len(leave_one_out_support.conflicting_clades)
    baseline_high_support_conflicts = _count_high_support_conflicts(baseline_support)
    leave_one_out_high_support_conflicts = _count_high_support_conflicts(
        leave_one_out_support
    )
    rooted_delta = (
        leave_one_out_topology.rooted_robinson_foulds_distance
        - baseline_topology.rooted_robinson_foulds_distance
    )
    unrooted_delta = (
        leave_one_out_topology.unrooted_robinson_foulds_distance
        - baseline_topology.unrooted_robinson_foulds_distance
    )
    rooted_normalized_delta = round(
        leave_one_out_topology.rooted_normalized_robinson_foulds
        - baseline_topology.rooted_normalized_robinson_foulds,
        6,
    )
    unrooted_normalized_delta = round(
        leave_one_out_topology.unrooted_normalized_robinson_foulds
        - baseline_topology.unrooted_normalized_robinson_foulds,
        6,
    )
    support_disagreement_delta = (
        leave_one_out_support_disagreements - baseline_support_disagreements
    )
    conflicting_clade_delta = (
        leave_one_out_conflicting_clades - baseline_conflicting_clades
    )
    high_support_conflict_delta = (
        leave_one_out_high_support_conflicts - baseline_high_support_conflicts
    )
    topology_changed = any(
        [
            leave_one_out_topology.topology_equal != baseline_topology.topology_equal,
            leave_one_out_topology.same_unrooted_topology
            != baseline_topology.same_unrooted_topology,
            leave_one_out_topology.same_taxa_different_rooting
            != baseline_topology.same_taxa_different_rooting,
            rooted_delta != 0,
            unrooted_delta != 0,
            rooted_normalized_delta != 0.0,
            unrooted_normalized_delta != 0.0,
        ]
    )
    support_changed = any(
        [
            support_disagreement_delta != 0,
            conflicting_clade_delta != 0,
            high_support_conflict_delta != 0,
        ]
    )
    influence_score = round(
        abs(rooted_normalized_delta)
        + abs(unrooted_normalized_delta)
        + abs(float(support_disagreement_delta))
        + abs(float(conflicting_clade_delta))
        + abs(float(high_support_conflict_delta)),
        6,
    )
    return TaxonInfluenceRow(
        taxon=taxon,
        retained_taxa=retained_taxa,
        baseline_topology_equal=baseline_topology.topology_equal,
        leave_one_out_topology_equal=leave_one_out_topology.topology_equal,
        baseline_same_unrooted_topology=baseline_topology.same_unrooted_topology,
        leave_one_out_same_unrooted_topology=(
            leave_one_out_topology.same_unrooted_topology
        ),
        baseline_same_taxa_different_rooting=(
            baseline_topology.same_taxa_different_rooting
        ),
        leave_one_out_same_taxa_different_rooting=(
            leave_one_out_topology.same_taxa_different_rooting
        ),
        baseline_rooted_robinson_foulds_distance=(
            baseline_topology.rooted_robinson_foulds_distance
        ),
        leave_one_out_rooted_robinson_foulds_distance=(
            leave_one_out_topology.rooted_robinson_foulds_distance
        ),
        rooted_robinson_foulds_delta=rooted_delta,
        baseline_rooted_normalized_robinson_foulds=(
            baseline_topology.rooted_normalized_robinson_foulds
        ),
        leave_one_out_rooted_normalized_robinson_foulds=(
            leave_one_out_topology.rooted_normalized_robinson_foulds
        ),
        rooted_normalized_robinson_foulds_delta=rooted_normalized_delta,
        baseline_unrooted_robinson_foulds_distance=(
            baseline_topology.unrooted_robinson_foulds_distance
        ),
        leave_one_out_unrooted_robinson_foulds_distance=(
            leave_one_out_topology.unrooted_robinson_foulds_distance
        ),
        unrooted_robinson_foulds_delta=unrooted_delta,
        baseline_unrooted_normalized_robinson_foulds=(
            baseline_topology.unrooted_normalized_robinson_foulds
        ),
        leave_one_out_unrooted_normalized_robinson_foulds=(
            leave_one_out_topology.unrooted_normalized_robinson_foulds
        ),
        unrooted_normalized_robinson_foulds_delta=unrooted_normalized_delta,
        baseline_support_disagreements=baseline_support_disagreements,
        leave_one_out_support_disagreements=leave_one_out_support_disagreements,
        support_disagreement_delta=support_disagreement_delta,
        baseline_conflicting_clades=baseline_conflicting_clades,
        leave_one_out_conflicting_clades=leave_one_out_conflicting_clades,
        conflicting_clade_delta=conflicting_clade_delta,
        baseline_high_support_conflicts=baseline_high_support_conflicts,
        leave_one_out_high_support_conflicts=leave_one_out_high_support_conflicts,
        high_support_conflict_delta=high_support_conflict_delta,
        topology_changed=topology_changed,
        support_changed=support_changed,
        influence_score=influence_score,
        influence_rank=0,
    )


def _rank_taxon_influence_rows(
    rows: list[TaxonInfluenceRow],
) -> list[TaxonInfluenceRow]:
    ranked = sorted(
        rows,
        key=lambda row: (
            -row.influence_score,
            not row.topology_changed,
            not row.support_changed,
            row.taxon,
        ),
    )
    for index, row in enumerate(ranked, start=1):
        row.influence_rank = index
    return ranked


def analyze_taxon_influence(
    left_path: Path,
    right_path: Path,
) -> TaxonInfluenceReport:
    """Measure how each shared taxon changes disagreement when excluded."""
    baseline_topology = _build_tree_comparison_report(
        left_path,
        right_path,
        _load_tree(left_path),
        _load_tree(right_path),
        rf_mode="rooted",
        taxon_overlap_policy="prune-to-shared",
    )
    if len(baseline_topology.shared_taxa) < 3:
        raise ValueError("taxon influence analysis requires at least three shared taxa")

    baseline_support = _build_support_comparison_report(
        left_path,
        right_path,
        _load_tree(left_path),
        _load_tree(right_path),
        strong_support_threshold=0.9,
        weak_support_threshold=0.7,
        support_disagreement_threshold=0.15,
    )
    shared_taxa = set(baseline_topology.shared_taxa)
    left_tree = _load_tree(left_path)
    right_tree = _load_tree(right_path)

    rows: list[TaxonInfluenceRow] = []
    for taxon in baseline_topology.shared_taxa:
        retained_taxa = sorted(shared_taxa - {taxon})
        retained_set = set(retained_taxa)
        left_pruned, _, _ = _prune_tree_against_shared_taxa(left_tree, retained_set)
        right_pruned, _, _ = _prune_tree_against_shared_taxa(right_tree, retained_set)
        leave_one_out_topology = _build_tree_comparison_report(
            left_path,
            right_path,
            left_pruned,
            right_pruned,
            rf_mode="rooted",
            taxon_overlap_policy="require-identical",
        )
        leave_one_out_support = _build_support_comparison_report(
            left_path,
            right_path,
            deepcopy(left_pruned),
            deepcopy(right_pruned),
            strong_support_threshold=0.9,
            weak_support_threshold=0.7,
            support_disagreement_threshold=0.15,
        )
        rows.append(
            _build_taxon_influence_row(
                taxon,
                retained_taxa,
                baseline_topology,
                leave_one_out_topology,
                baseline_support,
                leave_one_out_support,
            )
        )
    return TaxonInfluenceReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=baseline_topology.shared_taxa,
        baseline_topology=baseline_topology,
        baseline_support=baseline_support,
        rows=_rank_taxon_influence_rows(rows),
    )


def _prune_tree_against_shared_taxa(tree, keep_taxa: set[str]):
    from bijux_phylogenetics.phylo.pruning import _prune_tree_against_taxa

    return _prune_tree_against_taxa(tree, keep_taxa)


def write_taxon_influence_table(
    path: Path,
    left_path: Path,
    right_path: Path,
) -> Path:
    """Write one row per shared taxon from leave-one-out tree comparison."""
    report = analyze_taxon_influence(left_path, right_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "influence_rank",
                "taxon",
                "retained_taxa",
                "retained_taxon_count",
                "baseline_topology_equal",
                "leave_one_out_topology_equal",
                "baseline_same_unrooted_topology",
                "leave_one_out_same_unrooted_topology",
                "baseline_same_taxa_different_rooting",
                "leave_one_out_same_taxa_different_rooting",
                "baseline_rooted_robinson_foulds_distance",
                "leave_one_out_rooted_robinson_foulds_distance",
                "rooted_robinson_foulds_delta",
                "baseline_rooted_normalized_robinson_foulds",
                "leave_one_out_rooted_normalized_robinson_foulds",
                "rooted_normalized_robinson_foulds_delta",
                "baseline_unrooted_robinson_foulds_distance",
                "leave_one_out_unrooted_robinson_foulds_distance",
                "unrooted_robinson_foulds_delta",
                "baseline_unrooted_normalized_robinson_foulds",
                "leave_one_out_unrooted_normalized_robinson_foulds",
                "unrooted_normalized_robinson_foulds_delta",
                "baseline_support_disagreements",
                "leave_one_out_support_disagreements",
                "support_disagreement_delta",
                "baseline_conflicting_clades",
                "leave_one_out_conflicting_clades",
                "conflicting_clade_delta",
                "baseline_high_support_conflicts",
                "leave_one_out_high_support_conflicts",
                "high_support_conflict_delta",
                "topology_changed",
                "support_changed",
                "influence_score",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.rows:
            writer.writerow(
                {
                    "influence_rank": row.influence_rank,
                    "taxon": row.taxon,
                    "retained_taxa": "|".join(row.retained_taxa),
                    "retained_taxon_count": len(row.retained_taxa),
                    "baseline_topology_equal": str(row.baseline_topology_equal).lower(),
                    "leave_one_out_topology_equal": str(
                        row.leave_one_out_topology_equal
                    ).lower(),
                    "baseline_same_unrooted_topology": str(
                        row.baseline_same_unrooted_topology
                    ).lower(),
                    "leave_one_out_same_unrooted_topology": str(
                        row.leave_one_out_same_unrooted_topology
                    ).lower(),
                    "baseline_same_taxa_different_rooting": str(
                        row.baseline_same_taxa_different_rooting
                    ).lower(),
                    "leave_one_out_same_taxa_different_rooting": str(
                        row.leave_one_out_same_taxa_different_rooting
                    ).lower(),
                    "baseline_rooted_robinson_foulds_distance": (
                        row.baseline_rooted_robinson_foulds_distance
                    ),
                    "leave_one_out_rooted_robinson_foulds_distance": (
                        row.leave_one_out_rooted_robinson_foulds_distance
                    ),
                    "rooted_robinson_foulds_delta": row.rooted_robinson_foulds_delta,
                    "baseline_rooted_normalized_robinson_foulds": (
                        row.baseline_rooted_normalized_robinson_foulds
                    ),
                    "leave_one_out_rooted_normalized_robinson_foulds": (
                        row.leave_one_out_rooted_normalized_robinson_foulds
                    ),
                    "rooted_normalized_robinson_foulds_delta": (
                        row.rooted_normalized_robinson_foulds_delta
                    ),
                    "baseline_unrooted_robinson_foulds_distance": (
                        row.baseline_unrooted_robinson_foulds_distance
                    ),
                    "leave_one_out_unrooted_robinson_foulds_distance": (
                        row.leave_one_out_unrooted_robinson_foulds_distance
                    ),
                    "unrooted_robinson_foulds_delta": row.unrooted_robinson_foulds_delta,
                    "baseline_unrooted_normalized_robinson_foulds": (
                        row.baseline_unrooted_normalized_robinson_foulds
                    ),
                    "leave_one_out_unrooted_normalized_robinson_foulds": (
                        row.leave_one_out_unrooted_normalized_robinson_foulds
                    ),
                    "unrooted_normalized_robinson_foulds_delta": (
                        row.unrooted_normalized_robinson_foulds_delta
                    ),
                    "baseline_support_disagreements": (
                        row.baseline_support_disagreements
                    ),
                    "leave_one_out_support_disagreements": (
                        row.leave_one_out_support_disagreements
                    ),
                    "support_disagreement_delta": row.support_disagreement_delta,
                    "baseline_conflicting_clades": row.baseline_conflicting_clades,
                    "leave_one_out_conflicting_clades": (
                        row.leave_one_out_conflicting_clades
                    ),
                    "conflicting_clade_delta": row.conflicting_clade_delta,
                    "baseline_high_support_conflicts": (
                        row.baseline_high_support_conflicts
                    ),
                    "leave_one_out_high_support_conflicts": (
                        row.leave_one_out_high_support_conflicts
                    ),
                    "high_support_conflict_delta": row.high_support_conflict_delta,
                    "topology_changed": str(row.topology_changed).lower(),
                    "support_changed": str(row.support_changed).lower(),
                    "influence_score": row.influence_score,
                }
            )
    return path
