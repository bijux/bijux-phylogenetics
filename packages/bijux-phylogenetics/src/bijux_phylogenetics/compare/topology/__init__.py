"""Topology comparison workflows, support audits, and reviewer-facing ledgers."""

from __future__ import annotations

from importlib import import_module

_PUBLIC_SURFACES = (
    (
        "agreement",
        (
            "approximate_maximum_agreement_subtree",
            "prune_trees_to_agreement_subtree",
        ),
    ),
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
        "clade_ages",
        ("compare_clade_ages",),
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
        "deep_coalescence",
        ("count_deep_coalescences",),
    ),
    (
        "reconciliation",
        ("reconcile_duplication_loss_transfer",),
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
            "AgreementSubtreeCandidateRow",
            "AgreementSubtreePruningReport",
            "MaximumAgreementSubtreeApproximationReport",
            "MaximumAgreementSubtreeSearchRow",
            "BranchLengthComparisonReport",
            "BranchLengthPair",
            "BranchScoreComparisonReport",
            "BranchScoreSplit",
            "BranchScoreStatus",
            "CladeAgeComparisonRow",
            "CladeChangeReport",
            "DeepCoalescenceBranchRow",
            "DeepCoalescenceReport",
            "DeepCoalescenceTaxonMapRow",
            "DuplicationLossTransferAssociationRow",
            "DuplicationLossTransferEventRow",
            "DuplicationLossTransferReport",
            "DateAwareTreeComparisonReport",
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
        "structural_parity",
        (
            "StructuralTreeParityReport",
            "StructuralTreeSetParityReport",
            "compare_tree_sets_structurally",
            "compare_tree_structurally",
        ),
    ),
    (
        "tables",
        (
            "write_agreement_subtree_pruning_table",
            "write_agreement_subtree_removed_taxa_table",
            "write_agreement_subtree_search_table",
            "write_maximum_agreement_subtree_pruning_table",
            "write_maximum_agreement_subtree_removed_taxa_table",
            "write_maximum_agreement_subtree_search_table",
            "write_clade_overlap_table",
            "write_deep_coalescence_branch_table",
            "write_deep_coalescence_taxon_map_table",
            "write_duplication_loss_transfer_event_table",
            "write_duplication_loss_transfer_taxon_map_table",
            "write_date_aware_tree_comparison_table",
            "write_shared_taxa_pruning_table",
            "write_shared_taxa_removed_taxa_table",
            "write_support_comparison_table",
            "write_tree_comparison_table",
        ),
    ),
)

__all__ = [name for _, names in _PUBLIC_SURFACES for name in names]

_NAME_TO_MODULE = {
    name: module_name for module_name, names in _PUBLIC_SURFACES for name in names
}


def __getattr__(name: str):
    """Resolve topology exports lazily from their owning submodules."""
    try:
        module_name = _NAME_TO_MODULE[name]
    except KeyError as error:
        raise AttributeError(
            f"module {__name__!r} has no attribute {name!r}"
        ) from error
    value = getattr(import_module(f"{__name__}.{module_name}"), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
