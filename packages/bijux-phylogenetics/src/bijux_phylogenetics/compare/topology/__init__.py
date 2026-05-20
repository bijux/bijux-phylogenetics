from __future__ import annotations

from bijux_phylogenetics.core.clade_sets import (
    canonical_clade_id,
)
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
from .support import compare_support_values
from .support import _build_support_comparison_report
from .tables import (
    write_clade_overlap_table,
    write_shared_taxa_pruning_table,
    write_shared_taxa_removed_taxa_table,
    write_support_comparison_table,
    write_tree_comparison_table,
)


def _format_clade_set(clades: set[frozenset[str]]) -> list[str]:
    return sorted(canonical_clade_id(clade) for clade in clades)


def _split_id(signature: frozenset[str]) -> str:
    return canonical_clade_id(signature)
