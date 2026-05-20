"""Topology comparison workflows, support audits, and reviewer-facing ledgers."""

from __future__ import annotations

from importlib import import_module

_PUBLIC_SURFACES = (
    (
        "branch_lengths",
        (
            "_build_branch_score_report",
            "_compare_branch_lengths_for_trees",
            "compare_branch_lengths",
            "compare_branch_score_distance",
        ),
    ),
    (
        "comparison",
        (
            "_build_tree_comparison_report",
            "_compare_tree_objects",
            "_resolve_shared_taxa",
            "_validate_rf_mode",
            "_validate_taxon_overlap_policy",
            "compare_robinson_foulds",
            "compare_tree_paths",
        ),
    ),
    (
        "distance",
        (
            "TopologyDistanceReport",
            "TopologyDistanceSplitRow",
            "compare_topology_distance",
            "compare_topology_distance_trees",
            "write_topology_distance_split_table",
        ),
    ),
    (
        "models",
        (
            "BranchLengthComparisonReport",
            "BranchLengthPair",
            "BranchScoreComparisonReport",
            "BranchScoreSplit",
            "BranchScoreStatus",
            "CladeChangeReport",
            "CladeOverlapComparisonReport",
            "CladeOverlapObservation",
            "CladeOverlapRow",
            "CladeSetComparisonReport",
            "CladeSupportPair",
            "InMemoryBranchLengthComparison",
            "InMemoryTopologyComparison",
            "RobinsonFouldsComparisonReport",
            "RobinsonFouldsMode",
            "SharedTaxaPruningReport",
            "SupportComparisonReport",
            "SupportConflictRow",
            "TaxonOverlapPolicy",
            "TreeCladeOverlapSummary",
            "TreeComparisonReport",
        ),
    ),
    (
        "overlap",
        (
            "compare_clade_overlap",
            "compare_clade_sets",
            "detect_clade_changes",
            "prune_trees_to_shared_taxa",
        ),
    ),
    (
        "support",
        (
            "_build_support_comparison_report",
            "compare_support_values",
        ),
    ),
    (
        "tables",
        (
            "write_clade_overlap_table",
            "write_shared_taxa_pruning_table",
            "write_shared_taxa_removed_taxa_table",
            "write_support_comparison_table",
            "write_tree_comparison_table",
        ),
    ),
)

__all__ = [name for _, names in _PUBLIC_SURFACES for name in names]

_NAME_TO_MODULE = {
    name: module_name
    for module_name, names in _PUBLIC_SURFACES
    for name in names
}


def __getattr__(name: str):
    """Resolve topology exports lazily from their owning submodules."""
    try:
        module_name = _NAME_TO_MODULE[name]
    except KeyError as error:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from error
    value = getattr(import_module(f"{__name__}.{module_name}"), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
