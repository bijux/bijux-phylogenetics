"""Tree-set loading, consensus, support, and distance workflows."""

from .budgets import (
    build_tree_set_budget_report,
    build_tree_set_workflow_budget,
    enforce_tree_set_tree_budget,
)
from .clade_compatibility import (
    compute_clade_compatibility_graph,
    write_clade_compatibility_edge_table,
    write_clade_compatibility_graph_dot,
    write_clade_compatibility_node_table,
)
from .clade_support import (
    compute_clade_frequency_table,
    compute_reference_tree_clade_support,
    write_clade_frequency_table,
    write_reference_tree_clade_support_table,
)
from .consensus import (
    compute_consensus_tree,
    compute_consensus_tree_with_threshold,
    compute_strict_consensus_tree,
    write_consensus_tree,
)
from .contracts import (
    CladeCompatibilityEdgeRow,
    CladeCompatibilityGraphReport,
    CladeCompatibilityNodeRow,
    CladeFrequency,
    CladeFrequencyReport,
    ConsensusTreeReport,
    GeneTreeQuartetConcordanceReport,
    GeneTreeQuartetConcordanceRow,
    MajorityRuleExtendedAcceptedCladeRow,
    MajorityRuleExtendedConsensusReport,
    MajorityRuleExtendedRejectedCladeRow,
    QuartetPuzzlingAssemblyRow,
    QuartetPuzzlingReport,
    QuartetTopologyScoreRow,
    TreeDistanceMatrixReport,
    TreeDistancePair,
    TreeSetCladeSupportReport,
    TreeSetCladeSupportRow,
    TreeSetProcessingSummary,
    TreeSetQuartetSupportReport,
    TreeSetQuartetSupportRow,
    TreeSetRecord,
    TreeSetReport,
    TreeSetWorkflowBudget,
    TreeSetWorkflowBudgetReport,
)
from .distances import compute_tree_distance_matrix, write_tree_distance_matrix
from .extended_consensus import (
    compute_majority_rule_extended_consensus,
    write_majority_rule_extended_consensus_artifacts,
    write_majority_rule_extended_consensus_inclusion_table,
    write_majority_rule_extended_consensus_rejected_conflict_table,
)
from .inventory import load_tree_set
from .quartet_concordance import (
    compute_gene_tree_quartet_concordance_factors,
    write_gene_tree_quartet_concordance_table,
)
from .quartet_puzzling import (
    build_quartet_puzzling_consensus,
    write_quartet_puzzling_artifacts,
)
from .quartet_support import (
    compute_reference_tree_quartet_support,
    write_reference_tree_quartet_support_table,
)

__all__ = [
    "CladeFrequency",
    "CladeFrequencyReport",
    "CladeCompatibilityEdgeRow",
    "CladeCompatibilityGraphReport",
    "CladeCompatibilityNodeRow",
    "ConsensusTreeReport",
    "GeneTreeQuartetConcordanceReport",
    "GeneTreeQuartetConcordanceRow",
    "MajorityRuleExtendedAcceptedCladeRow",
    "MajorityRuleExtendedConsensusReport",
    "MajorityRuleExtendedRejectedCladeRow",
    "TreeDistanceMatrixReport",
    "TreeDistancePair",
    "TreeSetCladeSupportReport",
    "TreeSetCladeSupportRow",
    "QuartetPuzzlingAssemblyRow",
    "QuartetPuzzlingReport",
    "QuartetTopologyScoreRow",
    "TreeSetQuartetSupportReport",
    "TreeSetQuartetSupportRow",
    "TreeSetProcessingSummary",
    "TreeSetRecord",
    "TreeSetReport",
    "TreeSetWorkflowBudget",
    "TreeSetWorkflowBudgetReport",
    "build_quartet_puzzling_consensus",
    "build_tree_set_budget_report",
    "build_tree_set_workflow_budget",
    "compute_clade_compatibility_graph",
    "compute_clade_frequency_table",
    "compute_consensus_tree",
    "compute_gene_tree_quartet_concordance_factors",
    "compute_majority_rule_extended_consensus",
    "compute_consensus_tree_with_threshold",
    "compute_reference_tree_clade_support",
    "compute_reference_tree_quartet_support",
    "compute_strict_consensus_tree",
    "compute_tree_distance_matrix",
    "enforce_tree_set_tree_budget",
    "load_tree_set",
    "write_clade_compatibility_edge_table",
    "write_clade_compatibility_graph_dot",
    "write_clade_compatibility_node_table",
    "write_clade_frequency_table",
    "write_consensus_tree",
    "write_gene_tree_quartet_concordance_table",
    "write_majority_rule_extended_consensus_artifacts",
    "write_majority_rule_extended_consensus_inclusion_table",
    "write_majority_rule_extended_consensus_rejected_conflict_table",
    "write_quartet_puzzling_artifacts",
    "write_reference_tree_clade_support_table",
    "write_reference_tree_quartet_support_table",
    "write_tree_distance_matrix",
]
