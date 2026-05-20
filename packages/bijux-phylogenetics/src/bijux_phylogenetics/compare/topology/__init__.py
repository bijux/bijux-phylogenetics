from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.core.clade_sets import (
    canonical_clade_id,
    informative_rooted_clade_nodes,
    node_support_value,
    split_sort_key,
)
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.diagnostics.validation import _load_tree
from bijux_phylogenetics.io.iqtree_support import support_fraction
from .branch_lengths import (
    _build_branch_score_report,
    _compare_branch_lengths_for_trees,
    compare_branch_lengths,
    compare_branch_score_distance,
)
from .comparison import (
    _build_tree_comparison_report,
    _compare_tree_objects,
    _resolve_shared_taxa,
    _validate_rf_mode,
    _validate_taxon_overlap_policy,
    compare_robinson_foulds,
    compare_tree_paths,
)
from .models import (
    BranchLengthComparisonReport,
    BranchLengthPair,
    BranchScoreComparisonReport,
    BranchScoreSplit,
    BranchScoreStatus,
    CladeChangeReport,
    CladeOverlapComparisonReport,
    CladeOverlapObservation,
    CladeOverlapRow,
    CladeSetComparisonReport,
    CladeSupportPair,
    InMemoryBranchLengthComparison,
    InMemoryTopologyComparison,
    RobinsonFouldsComparisonReport,
    RobinsonFouldsMode,
    SharedTaxaPruningReport,
    SupportComparisonReport,
    SupportConflictRow,
    TaxonOverlapPolicy,
    TreeCladeOverlapSummary,
    TreeComparisonReport,
    _STRONG_SUPPORT_THRESHOLD,
    _SUPPORT_DISAGREEMENT_THRESHOLD,
    _WEAK_SUPPORT_THRESHOLD,
)
from .overlap import (
    compare_clade_overlap,
    compare_clade_sets,
    detect_clade_changes,
    prune_trees_to_shared_taxa,
)


def _format_clade_set(clades: set[frozenset[str]]) -> list[str]:
    return sorted(canonical_clade_id(clade) for clade in clades)


def _split_id(signature: frozenset[str]) -> str:
    return canonical_clade_id(signature)


def _support_strength(
    value: float | None,
    *,
    strong_support_threshold: float,
    weak_support_threshold: float,
) -> str:
    fraction = support_fraction(value)
    if fraction is None:
        return "unavailable"
    if fraction >= strong_support_threshold:
        return "strong"
    if fraction >= weak_support_threshold:
        return "moderate"
    return "low"


def _pipe_join(values: list[str]) -> str:
    return "|".join(values)


def write_shared_taxa_pruning_table(
    path: Path,
    left_path: Path,
    right_path: Path,
) -> Path:
    """Write one row per tree summarizing shared-taxon pruning evidence."""
    _, _, report = prune_trees_to_shared_taxa(left_path, right_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tree_side",
                "tree_path",
                "original_tip_count",
                "retained_tip_count",
                "removed_tip_count",
                "requested_taxa",
                "kept_taxa",
                "removed_taxa",
                "absent_requested_taxa",
                "removed_taxa_with_reasons",
                "transformation",
                "root_to_tip_complete",
                "min_root_to_tip",
                "max_root_to_tip",
                "unary_internal_nodes",
                "original_total_branch_length",
                "pruned_total_branch_length",
                "branch_length_delta",
                "lost_taxa_count",
                "lost_taxa_fraction",
                "lost_clade_count",
                "lost_clade_fraction",
                "lost_branch_length",
                "lost_branch_length_fraction",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for tree_side, pruning in (
            ("left", report.left_pruning),
            ("right", report.right_pruning),
        ):
            writer.writerow(
                {
                    "tree_side": tree_side,
                    "tree_path": str(pruning.tree_path),
                    "original_tip_count": pruning.original_tip_count,
                    "retained_tip_count": len(pruning.kept_taxa),
                    "removed_tip_count": len(pruning.removed_taxa),
                    "requested_taxa": _pipe_join(pruning.requested_taxa),
                    "kept_taxa": _pipe_join(pruning.kept_taxa),
                    "removed_taxa": _pipe_join(pruning.removed_taxa),
                    "absent_requested_taxa": _pipe_join(pruning.absent_requested_taxa),
                    "removed_taxa_with_reasons": _pipe_join(
                        [
                            f"{row.taxon}:{row.reason}"
                            for row in pruning.removed_taxa_with_reasons
                        ]
                    ),
                    "transformation": pruning.summary.transformation,
                    "root_to_tip_complete": str(
                        pruning.pruning_audit.root_to_tip_complete
                    ).lower(),
                    "min_root_to_tip": pruning.pruning_audit.min_root_to_tip,
                    "max_root_to_tip": pruning.pruning_audit.max_root_to_tip,
                    "unary_internal_nodes": _pipe_join(
                        pruning.pruning_audit.unary_internal_nodes
                    ),
                    "original_total_branch_length": (
                        pruning.pruning_audit.original_total_branch_length
                    ),
                    "pruned_total_branch_length": (
                        pruning.pruning_audit.pruned_total_branch_length
                    ),
                    "branch_length_delta": pruning.pruning_audit.branch_length_delta,
                    "lost_taxa_count": pruning.information_loss.lost_taxa_count,
                    "lost_taxa_fraction": pruning.information_loss.lost_taxa_fraction,
                    "lost_clade_count": pruning.information_loss.lost_clade_count,
                    "lost_clade_fraction": pruning.information_loss.lost_clade_fraction,
                    "lost_branch_length": pruning.information_loss.lost_branch_length,
                    "lost_branch_length_fraction": (
                        pruning.information_loss.lost_branch_length_fraction
                    ),
                }
            )
    return path


def write_shared_taxa_removed_taxa_table(
    path: Path,
    left_path: Path,
    right_path: Path,
) -> Path:
    """Write one row per removed taxon from shared-taxon pruning."""
    _, _, report = prune_trees_to_shared_taxa(left_path, right_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["tree_side", "tree_path", "taxon", "reason"],
            delimiter="\t",
        )
        writer.writeheader()
        for tree_side, pruning in (
            ("left", report.left_pruning),
            ("right", report.right_pruning),
        ):
            for removed in pruning.removed_taxa_with_reasons:
                writer.writerow(
                    {
                        "tree_side": tree_side,
                        "tree_path": str(pruning.tree_path),
                        "taxon": removed.taxon,
                        "reason": removed.reason,
                    }
                )
    return path


def compare_support_values(
    left_path: Path,
    right_path: Path,
    *,
    strong_support_threshold: float = _STRONG_SUPPORT_THRESHOLD,
    weak_support_threshold: float = _WEAK_SUPPORT_THRESHOLD,
    support_disagreement_threshold: float = _SUPPORT_DISAGREEMENT_THRESHOLD,
) -> SupportComparisonReport:
    """Compare clade support values and support-aware conflicts across two trees."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    return _build_support_comparison_report(
        left_path,
        right_path,
        left,
        right,
        strong_support_threshold=strong_support_threshold,
        weak_support_threshold=weak_support_threshold,
        support_disagreement_threshold=support_disagreement_threshold,
    )


def _build_support_comparison_report(
    left_path: Path,
    right_path: Path,
    left: PhyloTree,
    right: PhyloTree,
    *,
    strong_support_threshold: float,
    weak_support_threshold: float,
    support_disagreement_threshold: float,
) -> SupportComparisonReport:
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa = left_taxa & right_taxa
    if len(shared_taxa) < 2:
        raise ValueError("support comparison requires at least two shared taxa")

    left_clades = informative_rooted_clade_nodes(left, shared_taxa)
    right_clades = informative_rooted_clade_nodes(right, shared_taxa)
    shared_clade_ids = sorted(
        left_clades.keys() & right_clades.keys(), key=split_sort_key
    )
    all_clade_ids = sorted(
        left_clades.keys() | right_clades.keys(),
        key=split_sort_key,
    )
    shared_clades: list[CladeSupportPair] = []
    conflicting_clades: list[SupportConflictRow] = []
    for clade_id in shared_clade_ids:
        left_support = node_support_value(left_clades[clade_id])
        right_support = node_support_value(right_clades[clade_id])
        left_fraction = support_fraction(left_support)
        right_fraction = support_fraction(right_support)
        shared_clades.append(
            CladeSupportPair(
                split_id=_split_id(clade_id),
                left_support=left_support,
                right_support=right_support,
                left_support_fraction=left_fraction,
                right_support_fraction=right_fraction,
                support_fraction_delta=(
                    None
                    if left_fraction is None or right_fraction is None
                    else abs(left_fraction - right_fraction)
                ),
                support_disagreement=(
                    left_fraction is not None
                    and right_fraction is not None
                    and abs(left_fraction - right_fraction)
                    >= support_disagreement_threshold
                ),
            )
        )
    for clade_id in all_clade_ids:
        left_present = clade_id in left_clades
        right_present = clade_id in right_clades
        if left_present and right_present:
            continue
        left_support = (
            node_support_value(left_clades[clade_id]) if left_present else None
        )
        right_support = (
            node_support_value(right_clades[clade_id]) if right_present else None
        )
        left_fraction = support_fraction(left_support)
        right_fraction = support_fraction(right_support)
        strongest_support_fraction = (
            max(value for value in (left_fraction, right_fraction) if value is not None)
            if left_fraction is not None or right_fraction is not None
            else None
        )
        support_strength = _support_strength(
            left_support if left_present else right_support,
            strong_support_threshold=strong_support_threshold,
            weak_support_threshold=weak_support_threshold,
        )
        if strongest_support_fraction is None:
            conflict_classification = "support_unavailable"
            detail = "clade conflict could not be ranked because no branch support was available"
        elif strongest_support_fraction >= strong_support_threshold:
            conflict_classification = "high_support_conflict"
            detail = "conflicting clade carried strong branch support in the tree where it was present"
        elif strongest_support_fraction >= weak_support_threshold:
            conflict_classification = "moderate_support_disagreement"
            detail = "conflicting clade carried moderate branch support in the tree where it was present"
        else:
            conflict_classification = "low_support_disagreement"
            detail = "conflicting clade was only weakly supported in the tree where it was present"
        conflicting_clades.append(
            SupportConflictRow(
                split_id=_split_id(clade_id),
                comparison_status="left_only" if left_present else "right_only",
                left_present=left_present,
                right_present=right_present,
                left_support=left_support,
                right_support=right_support,
                left_support_fraction=left_fraction,
                right_support_fraction=right_fraction,
                strongest_support_fraction=strongest_support_fraction,
                support_strength=support_strength,
                conflict_classification=conflict_classification,
                detail=detail,
            )
        )
    return SupportComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=sorted(shared_taxa),
        strong_support_threshold=strong_support_threshold,
        weak_support_threshold=weak_support_threshold,
        support_disagreement_threshold=support_disagreement_threshold,
        shared_clades=shared_clades,
        conflicting_clades=conflicting_clades,
    )


def write_clade_overlap_table(path: Path, tree_paths: list[Path]) -> Path:
    """Write one row per clade-per-tree overlap observation."""
    report = compare_clade_overlap(tree_paths)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "clade_id",
                "tree_path",
                "present",
                "support",
                "present_in_all_trees",
                "present_tree_count",
                "absent_tree_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.clade_rows:
            for observation in row.observations:
                writer.writerow(
                    {
                        "clade_id": row.clade_id,
                        "tree_path": str(observation.tree_path),
                        "present": str(observation.present).lower(),
                        "support": ""
                        if observation.support is None
                        else observation.support,
                        "present_in_all_trees": str(row.present_in_all_trees).lower(),
                        "present_tree_count": row.present_tree_count,
                        "absent_tree_count": row.absent_tree_count,
                    }
                )
    return path


def write_support_comparison_table(
    path: Path, left_path: Path, right_path: Path
) -> Path:
    """Write a flat TSV ledger for support-aware clade comparison."""
    report = compare_support_values(left_path, right_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "split_id",
                "row_kind",
                "comparison_status",
                "left_present",
                "right_present",
                "left_support",
                "right_support",
                "left_support_fraction",
                "right_support_fraction",
                "support_fraction_delta",
                "support_disagreement",
                "strongest_support_fraction",
                "support_strength",
                "conflict_classification",
                "detail",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.shared_clades:
            writer.writerow(
                {
                    "split_id": row.split_id,
                    "row_kind": "shared_clade",
                    "comparison_status": "shared",
                    "left_present": "true",
                    "right_present": "true",
                    "left_support": ""
                    if row.left_support is None
                    else row.left_support,
                    "right_support": ""
                    if row.right_support is None
                    else row.right_support,
                    "left_support_fraction": (
                        ""
                        if row.left_support_fraction is None
                        else row.left_support_fraction
                    ),
                    "right_support_fraction": (
                        ""
                        if row.right_support_fraction is None
                        else row.right_support_fraction
                    ),
                    "support_fraction_delta": (
                        ""
                        if row.support_fraction_delta is None
                        else row.support_fraction_delta
                    ),
                    "support_disagreement": str(row.support_disagreement).lower(),
                    "strongest_support_fraction": "",
                    "support_strength": "",
                    "conflict_classification": "",
                    "detail": (
                        "shared clade support differs across trees"
                        if row.support_disagreement
                        else "shared clade support is aligned across trees"
                    ),
                }
            )
        for row in report.conflicting_clades:
            writer.writerow(
                {
                    "split_id": row.split_id,
                    "row_kind": "conflicting_clade",
                    "comparison_status": row.comparison_status,
                    "left_present": str(row.left_present).lower(),
                    "right_present": str(row.right_present).lower(),
                    "left_support": ""
                    if row.left_support is None
                    else row.left_support,
                    "right_support": ""
                    if row.right_support is None
                    else row.right_support,
                    "left_support_fraction": (
                        ""
                        if row.left_support_fraction is None
                        else row.left_support_fraction
                    ),
                    "right_support_fraction": (
                        ""
                        if row.right_support_fraction is None
                        else row.right_support_fraction
                    ),
                    "support_fraction_delta": "",
                    "support_disagreement": "false",
                    "strongest_support_fraction": (
                        ""
                        if row.strongest_support_fraction is None
                        else row.strongest_support_fraction
                    ),
                    "support_strength": row.support_strength,
                    "conflict_classification": row.conflict_classification,
                    "detail": row.detail,
                }
            )
    return path


def write_tree_comparison_table(path: Path, left_path: Path, right_path: Path) -> Path:
    """Write a flat TSV table covering the compared clade and split surfaces."""
    clades = compare_clade_sets(left_path, right_path)
    support = compare_support_values(left_path, right_path)
    branch_lengths = compare_branch_lengths(left_path, right_path)
    support_by_id = {row.split_id: row for row in support.shared_clades}
    support_conflict_by_id = {row.split_id: row for row in support.conflicting_clades}
    branch_by_id = {row.split_id: row for row in branch_lengths.shared_splits}
    branch_score_by_id = {
        row.split_id: row for row in branch_lengths.branch_score.splits
    }
    all_split_ids = sorted(
        set(clades.shared_clades)
        | set(clades.left_only_clades)
        | set(clades.right_only_clades)
        | set(support_by_id)
        | set(support_conflict_by_id)
        | set(branch_by_id)
        | set(branch_score_by_id)
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "split_id",
                "comparison_status",
                "shared_clade",
                "left_support",
                "right_support",
                "left_support_fraction",
                "right_support_fraction",
                "support_fraction_delta",
                "support_disagreement",
                "support_conflict_classification",
                "support_conflict_strength",
                "left_length",
                "right_length",
                "length_delta",
                "length_ratio",
                "branch_score_status",
                "branch_score_difference",
                "branch_score_squared_difference",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for split_id in all_split_ids:
            support_row = support_by_id.get(split_id)
            support_conflict_row = support_conflict_by_id.get(split_id)
            branch_row = branch_by_id.get(split_id)
            branch_score_row = branch_score_by_id.get(split_id)
            if split_id in clades.shared_clades:
                status = "shared"
            elif split_id in clades.left_only_clades:
                status = "left_only"
            else:
                status = "right_only"
            writer.writerow(
                {
                    "split_id": split_id,
                    "comparison_status": status,
                    "shared_clade": str(split_id in clades.shared_clades).lower(),
                    "left_support": ""
                    if (
                        (support_row is None or support_row.left_support is None)
                        and (
                            support_conflict_row is None
                            or support_conflict_row.left_support is None
                        )
                    )
                    else (
                        support_row.left_support
                        if support_row is not None
                        else support_conflict_row.left_support
                    ),
                    "right_support": ""
                    if (
                        (support_row is None or support_row.right_support is None)
                        and (
                            support_conflict_row is None
                            or support_conflict_row.right_support is None
                        )
                    )
                    else (
                        support_row.right_support
                        if support_row is not None
                        else support_conflict_row.right_support
                    ),
                    "left_support_fraction": ""
                    if (
                        (
                            support_row is None
                            or support_row.left_support_fraction is None
                        )
                        and (
                            support_conflict_row is None
                            or support_conflict_row.left_support_fraction is None
                        )
                    )
                    else (
                        support_row.left_support_fraction
                        if support_row is not None
                        else support_conflict_row.left_support_fraction
                    ),
                    "right_support_fraction": ""
                    if (
                        (
                            support_row is None
                            or support_row.right_support_fraction is None
                        )
                        and (
                            support_conflict_row is None
                            or support_conflict_row.right_support_fraction is None
                        )
                    )
                    else (
                        support_row.right_support_fraction
                        if support_row is not None
                        else support_conflict_row.right_support_fraction
                    ),
                    "support_fraction_delta": ""
                    if support_row is None or support_row.support_fraction_delta is None
                    else support_row.support_fraction_delta,
                    "support_disagreement": (
                        "false"
                        if support_row is None
                        else str(support_row.support_disagreement).lower()
                    ),
                    "support_conflict_classification": (
                        ""
                        if support_conflict_row is None
                        else support_conflict_row.conflict_classification
                    ),
                    "support_conflict_strength": (
                        ""
                        if support_conflict_row is None
                        else support_conflict_row.support_strength
                    ),
                    "left_length": ""
                    if branch_row is None or branch_row.left_length is None
                    else branch_row.left_length,
                    "right_length": ""
                    if branch_row is None or branch_row.right_length is None
                    else branch_row.right_length,
                    "length_delta": ""
                    if branch_row is None or branch_row.delta is None
                    else branch_row.delta,
                    "length_ratio": ""
                    if branch_row is None or branch_row.ratio is None
                    else branch_row.ratio,
                    "branch_score_status": ""
                    if branch_score_row is None
                    else branch_score_row.comparison_status,
                    "branch_score_difference": ""
                    if branch_score_row is None
                    or branch_score_row.branch_score_difference is None
                    else branch_score_row.branch_score_difference,
                    "branch_score_squared_difference": ""
                    if branch_score_row is None
                    or branch_score_row.squared_difference is None
                    else branch_score_row.squared_difference,
                }
            )
    return path
