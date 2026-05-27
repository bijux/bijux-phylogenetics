from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS = {
    "RobinsonFouldsMetrics": "clades",
    "canonical_bipartition": "clades",
    "canonical_clade_id": "clades",
    "informative_rooted_clade_nodes": "clades",
    "informative_rooted_clades": "clades",
    "informative_unrooted_splits": "clades",
    "node_support_value": "clades",
    "robinson_foulds_metrics": "clades",
    "split_sort_key": "clades",
    "tree_has_polytomy": "clades",
    "ape_node_id_for_node": "node_identity",
    "build_ape_internal_node_map": "node_identity",
    "build_ape_tip_node_map": "node_identity",
    "iter_internal_nodes_preorder": "node_identity",
    "BranchCollapseReport": "models",
    "CladeExtractionReport": "models",
    "SubtreeExtractionReport": "models",
    "TreeMonophylyReport": "models",
    "TreeMrcaReport": "models",
    "TreeOrderingReport": "models",
    "TreeRootingReport": "models",
    "TreeTransformationSummary": "models",
    "build_bionj_tree": "bionj",
    "build_neighbor_joining_tree": "neighbor_joining",
    "_order_tree": "ordering",
    "_rotate_all_nodes": "ordering",
    "_rotate_named_node": "ordering",
    "ladderize_tree": "ordering",
    "rotate_all_internal_nodes": "ordering",
    "rotate_named_node": "ordering",
    "sort_tree_tips_alphabetically": "ordering",
    "_adjacent_nodes": "rooting",
    "_analyze_midpoint_path": "rooting",
    "_biophylo_clade_taxa": "rooting",
    "_clone_subtree_away_from": "rooting",
    "_clone_subtree_component": "rooting",
    "_copy_node_payload": "rooting",
    "_edge_length_between": "rooting",
    "_find_monophyletic_outgroup_node": "rooting",
    "_is_strictly_bifurcating": "rooting",
    "_normalize_outgroup_rooting_to_ape": "rooting",
    "_root_tree_by_outgroup_node": "rooting",
    "reroot_tree_by_midpoint": "rooting",
    "root_tree_on_outgroup": "rooting",
    "tree_from_biophylo": "rooting",
    "tree_to_biophylo": "rooting",
    "unroot_tree": "rooting",
    "write_tree_rooting_report": "rooting",
    "_build_subtree": "subtree",
    "_common_prefix_length": "subtree",
    "_extract_subtree_report": "subtree",
    "_find_named_nodes": "subtree",
    "_interpreted_rooted_state": "subtree",
    "_monophyly_report_from_node": "subtree",
    "_mrca_node_from_taxa": "subtree",
    "_root_to_tip_paths": "subtree",
    "assess_tree_monophyly": "subtree",
    "extract_named_clade": "subtree",
    "extract_tree_clade_by_descendant_taxa": "subtree",
    "extract_tree_clade_by_node_id": "subtree",
    "find_tree_mrca": "subtree",
    "TipDistanceMatrixReport": "tip_distances",
    "TipDistanceMatrixRow": "tip_distances",
    "compute_tree_tip_distance_matrix": "tip_distances",
    "summarize_tree_tip_distances": "tip_distances",
    "write_tree_tip_distance_long_table": "tip_distances",
    "write_tree_tip_distance_matrix": "tip_distances",
    "_clone_node": "transformation",
    "_collapse_short_internal_branches": "transformation",
    "_combine_branch_lengths": "transformation",
    "_compare_tree_topology": "transformation",
    "_copy_node_without_children": "transformation",
    "_descendant_taxa": "transformation",
    "_format_optional_float": "transformation",
    "_join_taxa": "transformation",
    "_leaf_count": "transformation",
    "_node_signature": "transformation",
    "_summarize_transformation": "transformation",
    "collapse_branches_below_length": "transformation",
    "PhyloTree": "tree",
    "TaxonLabel": "tree",
    "TreeNode": "tree",
    "descendant_taxa": "tree",
    "normalize_taxon_key": "tree",
    "stable_node_label": "tree",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str) -> Any:
    module_name = _EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(f".{module_name}", __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
