from __future__ import annotations

import csv
from pathlib import Path

from .branch_lengths import compare_branch_lengths
from .agreement import (
    approximate_maximum_agreement_subtree,
    prune_trees_to_agreement_subtree,
)
from .overlap import (
    compare_clade_overlap,
    compare_clade_sets,
    prune_trees_to_shared_taxa,
)
from .models import RobinsonFouldsMode
from .support import compare_support_values


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


def write_agreement_subtree_pruning_table(
    path: Path,
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
) -> Path:
    """Write one row per tree summarizing agreement-subtree pruning evidence."""
    _, _, report = prune_trees_to_agreement_subtree(
        left_path,
        right_path,
        rf_mode=rf_mode,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tree_side",
                "tree_path",
                "rf_mode",
                "search_strategy",
                "possible_retained_subset_count",
                "evaluated_candidate_count",
                "agreement_retained_taxa",
                "agreement_removed_taxa",
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
                    "rf_mode": report.rf_mode,
                    "search_strategy": report.search_strategy,
                    "possible_retained_subset_count": report.possible_retained_subset_count,
                    "evaluated_candidate_count": report.evaluated_candidate_count,
                    "agreement_retained_taxa": _pipe_join(report.retained_taxa),
                    "agreement_removed_taxa": _pipe_join(report.agreement_removed_taxa),
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


def write_agreement_subtree_removed_taxa_table(
    path: Path,
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
) -> Path:
    """Write one row per taxon removed by agreement-subtree pruning."""
    _, _, report = prune_trees_to_agreement_subtree(
        left_path,
        right_path,
        rf_mode=rf_mode,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tree_side",
                "tree_path",
                "taxon",
                "reason",
                "shared_taxon",
                "removed_for_agreement_subtree",
            ],
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
                        "shared_taxon": str(
                            removed.taxon in report.shared_taxa
                        ).lower(),
                        "removed_for_agreement_subtree": str(
                            removed.taxon in report.agreement_removed_taxa
                        ).lower(),
                    }
                )
    return path


def write_agreement_subtree_search_table(
    path: Path,
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
) -> Path:
    """Write one row per evaluated retained-taxon candidate subset."""
    _, _, report = prune_trees_to_agreement_subtree(
        left_path,
        right_path,
        rf_mode=rf_mode,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "candidate_index",
                "retained_taxon_count",
                "retained_taxa",
                "removed_taxa",
                "robinson_foulds_distance",
                "normalized_robinson_foulds",
                "topology_equal",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.candidate_rows:
            writer.writerow(
                {
                    "candidate_index": row.candidate_index,
                    "retained_taxon_count": row.retained_taxon_count,
                    "retained_taxa": _pipe_join(row.retained_taxa),
                    "removed_taxa": _pipe_join(row.removed_taxa),
                    "robinson_foulds_distance": row.robinson_foulds_distance,
                    "normalized_robinson_foulds": row.normalized_robinson_foulds,
                    "topology_equal": str(row.topology_equal).lower(),
                }
            )
    return path


def write_maximum_agreement_subtree_pruning_table(
    path: Path,
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
    max_evaluated_candidate_count: int,
) -> Path:
    """Write one row per tree summarizing heuristic maximum-agreement pruning."""
    _, _, report = approximate_maximum_agreement_subtree(
        left_path,
        right_path,
        rf_mode=rf_mode,
        max_evaluated_candidate_count=max_evaluated_candidate_count,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tree_side",
                "tree_path",
                "rf_mode",
                "search_strategy",
                "selection_objective",
                "approximation_status",
                "possible_retained_subset_count",
                "max_evaluated_candidate_count",
                "evaluated_candidate_count",
                "maximum_agreement_retained_taxa",
                "maximum_agreement_removed_taxa",
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
                    "rf_mode": report.rf_mode,
                    "search_strategy": report.search_strategy,
                    "selection_objective": report.selection_objective,
                    "approximation_status": report.approximation_status,
                    "possible_retained_subset_count": report.possible_retained_subset_count,
                    "max_evaluated_candidate_count": report.max_evaluated_candidate_count,
                    "evaluated_candidate_count": report.evaluated_candidate_count,
                    "maximum_agreement_retained_taxa": _pipe_join(report.retained_taxa),
                    "maximum_agreement_removed_taxa": _pipe_join(
                        report.approximation_removed_taxa
                    ),
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


def write_maximum_agreement_subtree_removed_taxa_table(
    path: Path,
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
    max_evaluated_candidate_count: int,
) -> Path:
    """Write one row per taxon removed by heuristic maximum-agreement pruning."""
    _, _, report = approximate_maximum_agreement_subtree(
        left_path,
        right_path,
        rf_mode=rf_mode,
        max_evaluated_candidate_count=max_evaluated_candidate_count,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tree_side",
                "tree_path",
                "taxon",
                "reason",
                "shared_taxon",
                "removed_for_maximum_agreement_subtree",
            ],
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
                        "shared_taxon": str(removed.taxon in report.shared_taxa).lower(),
                        "removed_for_maximum_agreement_subtree": str(
                            removed.taxon in report.approximation_removed_taxa
                        ).lower(),
                    }
                )
    return path


def write_maximum_agreement_subtree_search_table(
    path: Path,
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
    max_evaluated_candidate_count: int,
) -> Path:
    """Write one row per heuristic candidate evaluated during maximum-agreement search."""
    _, _, report = approximate_maximum_agreement_subtree(
        left_path,
        right_path,
        rf_mode=rf_mode,
        max_evaluated_candidate_count=max_evaluated_candidate_count,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "evaluation_index",
                "step_index",
                "retained_taxon_count",
                "retained_taxa",
                "removed_taxa",
                "robinson_foulds_distance",
                "normalized_robinson_foulds",
                "topology_equal",
                "selected_for_next_step",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.search_rows:
            writer.writerow(
                {
                    "evaluation_index": row.evaluation_index,
                    "step_index": row.step_index,
                    "retained_taxon_count": row.retained_taxon_count,
                    "retained_taxa": _pipe_join(row.retained_taxa),
                    "removed_taxa": _pipe_join(row.removed_taxa),
                    "robinson_foulds_distance": row.robinson_foulds_distance,
                    "normalized_robinson_foulds": row.normalized_robinson_foulds,
                    "topology_equal": str(row.topology_equal).lower(),
                    "selected_for_next_step": str(row.selected_for_next_step).lower(),
                }
            )
    return path


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
