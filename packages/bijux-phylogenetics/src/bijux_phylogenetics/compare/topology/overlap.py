from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.diagnostics.validation.structure import _load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.clades import (
    informative_rooted_clade_nodes,
    informative_rooted_clades,
    node_support_value,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .comparison import _build_tree_comparison_report
from .models import (
    CladeChangeReport,
    CladeOverlapComparisonReport,
    CladeOverlapObservation,
    CladeOverlapRow,
    CladeSetComparisonReport,
    SharedTaxaPruningReport,
    TreeCladeOverlapSummary,
)


def _split_id(signature: frozenset[str]) -> str:
    return "|".join(sorted(signature))


def _resolve_shared_taxa_for_many_trees(
    tree_paths: list[Path],
) -> tuple[list[PhyloTree], set[str], list[list[str]]]:
    if len(tree_paths) < 2:
        raise ValueError("clade-overlap comparison requires at least two trees")
    trees = [_load_tree(path) for path in tree_paths]
    taxon_sets = [set(tree.tip_names) for tree in trees]
    shared_taxa = set.intersection(*taxon_sets)
    if len(shared_taxa) < 2:
        raise ValueError("clade-overlap comparison requires at least two shared taxa")
    excluded_taxa = [sorted(taxa - shared_taxa) for taxa in taxon_sets]
    return trees, shared_taxa, excluded_taxa


def prune_trees_to_shared_taxa(
    left_path: Path,
    right_path: Path,
) -> tuple[PhyloTree, PhyloTree, SharedTaxaPruningReport]:
    """Prune two trees to the exact shared taxon set."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa = sorted(left_taxa & right_taxa)
    if len(shared_taxa) < 2:
        raise ValueError("shared-taxon pruning requires at least two shared taxa")

    pruned_left, left_pruning = prune_tree_to_requested_taxa(left_path, shared_taxa)
    pruned_right, right_pruning = prune_tree_to_requested_taxa(right_path, shared_taxa)
    post_pruning_comparison = _build_tree_comparison_report(
        left_path,
        right_path,
        pruned_left,
        pruned_right,
        rf_mode="rooted",
        taxon_overlap_policy="require-identical",
    )
    return (
        pruned_left,
        pruned_right,
        SharedTaxaPruningReport(
            left_path=left_path,
            right_path=right_path,
            shared_taxa=shared_taxa,
            left_only_taxa=sorted(left_taxa - right_taxa),
            right_only_taxa=sorted(right_taxa - left_taxa),
            left_pruning=left_pruning,
            right_pruning=right_pruning,
            post_pruning_comparison=post_pruning_comparison,
        ),
    )


def compare_clade_sets(left_path: Path, right_path: Path) -> CladeSetComparisonReport:
    """Compare rooted informative clade sets across two trees."""
    overlap = compare_clade_overlap([left_path, right_path])
    left_summary, right_summary = overlap.tree_summaries
    return CladeSetComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=overlap.shared_taxa,
        shared_clades=overlap.shared_clades,
        left_only_clades=left_summary.unique_clades,
        right_only_clades=right_summary.unique_clades,
    )


def compare_clade_overlap(tree_paths: list[Path]) -> CladeOverlapComparisonReport:
    """Compare rooted clade overlap across two or more trees."""
    trees, shared_taxa, excluded_taxa = _resolve_shared_taxa_for_many_trees(tree_paths)
    clade_maps = [informative_rooted_clades(tree, shared_taxa) for tree in trees]
    clade_node_maps = [
        informative_rooted_clade_nodes(tree, shared_taxa) for tree in trees
    ]
    all_clades = sorted(
        set().union(*clade_maps),
        key=lambda signature: (len(signature), tuple(sorted(signature))),
    )
    shared_clades = [
        _split_id(clade)
        for clade in all_clades
        if all(clade in clade_map for clade_map in clade_maps)
    ]
    conflicting_clades = [
        _split_id(clade)
        for clade in all_clades
        if not all(clade in clade_map for clade_map in clade_maps)
    ]
    tree_summaries: list[TreeCladeOverlapSummary] = []
    for path, clade_map, clade_node_map, tree_excluded_taxa in zip(
        tree_paths,
        clade_maps,
        clade_node_maps,
        excluded_taxa,
        strict=True,
    ):
        unique_clades = [
            _split_id(clade)
            for clade in sorted(
                (
                    clade
                    for clade in clade_map
                    if sum(clade in other_map for other_map in clade_maps) == 1
                ),
                key=lambda signature: (len(signature), tuple(sorted(signature))),
            )
        ]
        support_clade_count = sum(
            1
            for clade in clade_map
            if node_support_value(clade_node_map[clade]) is not None
        )
        tree_summaries.append(
            TreeCladeOverlapSummary(
                tree_path=path,
                clade_count=len(clade_map),
                support_clade_count=support_clade_count,
                unique_clades=unique_clades,
                excluded_taxa=tree_excluded_taxa,
            )
        )
    clade_rows: list[CladeOverlapRow] = []
    for clade in all_clades:
        observations: list[CladeOverlapObservation] = []
        present_tree_count = 0
        for path, clade_map, clade_node_map in zip(
            tree_paths, clade_maps, clade_node_maps, strict=True
        ):
            present = clade in clade_map
            if present:
                present_tree_count += 1
            support = None
            if present:
                support = node_support_value(clade_node_map[clade])
            observations.append(
                CladeOverlapObservation(
                    tree_path=path,
                    present=present,
                    support=support,
                )
            )
        clade_rows.append(
            CladeOverlapRow(
                clade_id=_split_id(clade),
                present_in_all_trees=present_tree_count == len(tree_paths),
                present_tree_count=present_tree_count,
                absent_tree_count=len(tree_paths) - present_tree_count,
                observations=observations,
            )
        )
    return CladeOverlapComparisonReport(
        tree_paths=tree_paths,
        shared_taxa=sorted(shared_taxa),
        shared_clades=shared_clades,
        conflicting_clades=conflicting_clades,
        tree_summaries=tree_summaries,
        clade_rows=clade_rows,
    )


def detect_clade_changes(left_path: Path, right_path: Path) -> CladeChangeReport:
    """Report clades lost from the left tree and gained in the right tree."""
    report = compare_clade_sets(left_path, right_path)
    return CladeChangeReport(
        left_path=left_path,
        right_path=right_path,
        lost_clades=report.left_only_clades,
        gained_clades=report.right_only_clades,
    )
