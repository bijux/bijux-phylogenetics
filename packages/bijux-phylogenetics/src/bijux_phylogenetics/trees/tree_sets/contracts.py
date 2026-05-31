from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class TreeSetRecord:
    index: int
    tip_count: int
    taxa: list[str]
    rooted_topology_id: str
    unrooted_topology_id: str


@dataclass(frozen=True, slots=True)
class TreeSetProcessingSummary:
    runtime_seconds: float
    peak_memory_bytes: int
    skipped_malformed_tree_count: int


@dataclass(frozen=True, slots=True)
class TreeSetWorkflowBudget:
    max_tree_count: int | None = None
    max_report_table_rows: int | None = None
    memory_warning_threshold_bytes: int | None = None


@dataclass(slots=True)
class TreeSetWorkflowBudgetReport:
    max_tree_count: int | None
    max_report_table_rows: int | None
    memory_warning_threshold_bytes: int | None
    truncated_section_names: list[str]
    warning_messages: list[str]


@dataclass(slots=True)
class TreeSetReport:
    path: Path
    source_format: str
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    taxa_union: list[str]
    rooted_topology_count: int
    unrooted_topology_count: int
    records: list[TreeSetRecord]


@dataclass(frozen=True, slots=True)
class CladeFrequency:
    clade: str
    tree_count: int
    frequency: float


@dataclass(slots=True)
class CladeFrequencyReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    clade_frequencies: list[CladeFrequency]


@dataclass(frozen=True, slots=True)
class TreeSetSplitFrequencyRow:
    split: str
    tree_count: int
    frequency: float


@dataclass(slots=True)
class TreeSetSplitFrequencyReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    split_policy: str
    split_frequencies: list[TreeSetSplitFrequencyRow]


@dataclass(frozen=True, slots=True)
class TreeSetMaximumCladeCredibilityCandidateRow:
    score_rank: int
    source_tree_index: int
    rooted_topology_id: str
    raw_tree_count: int
    raw_tree_frequency: float
    clade_credibility_score: float
    candidate_newick: str


@dataclass(slots=True)
class TreeSetMaximumCladeCredibilityReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    rooted_topology_count: int
    selected_tree_index: int
    selected_rooted_topology_id: str
    clade_credibility_score: float
    maximum_clade_credibility_newick: str
    rows: list[TreeSetMaximumCladeCredibilityCandidateRow]


@dataclass(frozen=True, slots=True)
class TreeSetCredibleCladeRow:
    inclusion_rank: int
    clade: str
    tree_count: int
    frequency: float
    cumulative_frequency: float


@dataclass(slots=True)
class TreeSetCredibleCladeSetReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    credible_threshold: float
    included_clade_count: int
    excluded_clade_count: int
    included_cumulative_frequency: float
    included_clades: list[TreeSetCredibleCladeRow]
    excluded_clades: list[TreeSetCredibleCladeRow]


@dataclass(frozen=True, slots=True)
class PosteriorBranchLengthSummaryRow:
    clade: str
    matched_tree_count: int
    posterior_tree_count: int
    clade_frequency: float
    mean_branch_length: float
    median_branch_length: float
    hpd_95_lower: float
    hpd_95_upper: float
    effective_sample_size: float | None


@dataclass(slots=True)
class PosteriorBranchLengthSummaryReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    hpd_mass: float
    rows: list[PosteriorBranchLengthSummaryRow]


@dataclass(frozen=True, slots=True)
class PosteriorNodeAgeSummaryRow:
    clade: str
    node_kind: str
    matched_tree_count: int
    posterior_tree_count: int
    clade_frequency: float
    mean_node_age: float
    median_node_age: float
    hpd_95_lower: float
    hpd_95_upper: float
    effective_sample_size: float | None


@dataclass(slots=True)
class PosteriorNodeAgeSummaryReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    hpd_mass: float
    rows: list[PosteriorNodeAgeSummaryRow]


@dataclass(frozen=True, slots=True)
class PosteriorAgreementSubtreeCandidateRow:
    candidate_index: int
    retained_taxon_count: int
    retained_taxa: list[str]
    removed_taxa: list[str]
    rooted_topology_count: int
    dominant_rooted_topology_frequency: float
    stable_topology_reached: bool


@dataclass(slots=True)
class PosteriorAgreementSubtreeReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    search_strategy: str
    possible_retained_subset_count: int
    evaluated_candidate_count: int
    retained_taxa: list[str]
    agreement_removed_taxa: list[str]
    stable_rooted_topology_id: str
    agreement_subtree_newick: str
    candidate_rows: list[PosteriorAgreementSubtreeCandidateRow]


@dataclass(frozen=True, slots=True)
class PosteriorCladeCorrelationRow:
    left_clade: str
    right_clade: str
    compatibility_relation: str
    compatibility_reason: str
    left_tree_count: int
    right_tree_count: int
    left_frequency: float
    right_frequency: float
    cooccurrence_tree_count: int
    cooccurrence_frequency: float
    expected_cooccurrence_frequency: float
    binary_correlation: float | None


@dataclass(slots=True)
class PosteriorCladeCorrelationReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    clade_count: int
    pair_count: int
    clade_order: list[str]
    rows: list[PosteriorCladeCorrelationRow]


@dataclass(frozen=True, slots=True)
class PosteriorTreeDistanceDiagnosticRow:
    source_tree_index: int
    rooted_topology_id: str
    mcc_outlier_rank: int
    consensus_outlier_rank: int
    mcc_robinson_foulds_distance: int
    mcc_normalized_robinson_foulds: float
    mcc_branch_score_distance: float | None
    consensus_robinson_foulds_distance: int
    consensus_normalized_robinson_foulds: float
    consensus_branch_score_distance: float | None


@dataclass(frozen=True, slots=True)
class PosteriorTreeDistanceDistributionRow:
    reference_tree_kind: str
    distance_metric: str
    observed_value: float
    tree_count: int
    frequency: float


@dataclass(slots=True)
class PosteriorTreeDistanceDiagnosticsReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    maximum_clade_credibility_tree_index: int
    maximum_clade_credibility_rooted_topology_id: str
    maximum_clade_credibility_newick: str
    consensus_method: str
    consensus_newick: str
    row_count: int
    distribution_row_count: int
    rows: list[PosteriorTreeDistanceDiagnosticRow]
    distributions: list[PosteriorTreeDistanceDistributionRow]


@dataclass(frozen=True, slots=True)
class CladeCompatibilityNodeRow:
    clade: str
    tree_count: int
    frequency: float
    compatible_neighbor_count: int
    conflict_neighbor_count: int


@dataclass(frozen=True, slots=True)
class CladeCompatibilityEdgeRow:
    left_clade: str
    right_clade: str
    compatibility_relation: str
    compatibility_reason: str
    left_tree_count: int
    right_tree_count: int
    left_frequency: float
    right_frequency: float


@dataclass(slots=True)
class CladeCompatibilityGraphReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    node_count: int
    edge_count: int
    compatible_edge_count: int
    conflict_edge_count: int
    nodes: list[CladeCompatibilityNodeRow]
    edges: list[CladeCompatibilityEdgeRow]


@dataclass(frozen=True, slots=True)
class MajorityRuleExtendedAcceptedCladeRow:
    insertion_rank: int
    clade: str
    tree_count: int
    frequency: float
    inclusion_stage: str


@dataclass(frozen=True, slots=True)
class MajorityRuleExtendedRejectedCladeRow:
    clade: str
    tree_count: int
    frequency: float
    blocking_clades: list[str]


@dataclass(slots=True)
class MajorityRuleExtendedConsensusReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    consensus_method: str
    majority_threshold: float
    included_clade_count: int
    majority_included_clade_count: int
    extension_included_clade_count: int
    rejected_conflict_count: int
    consensus_newick: str
    accepted_clades: list[MajorityRuleExtendedAcceptedCladeRow]
    rejected_clades: list[MajorityRuleExtendedRejectedCladeRow]


@dataclass(frozen=True, slots=True)
class TreeSetCladeSupportRow:
    node_id: int
    node_kind: str
    node_label: str | None
    descendant_taxa: list[str]
    supporting_tree_count: int | None
    clade_frequency: float | None
    support_percent: float | None
    support_status: str
    explanation: str
    reference_branch_length: float | None
    reference_root_depth: float | None


@dataclass(slots=True)
class TreeSetCladeSupportReport:
    reference_tree_path: Path
    comparison_tree_set_path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    supported_clade_count: int
    absent_clade_count: int
    unscored_clade_count: int
    rows: list[TreeSetCladeSupportRow]


@dataclass(slots=True)
class TreeSetQuartetSupportRow:
    branch_id: str
    left_taxa: list[str]
    right_taxa: list[str]
    quartet_count_per_tree: int
    concordant_quartet_count: int
    discordant_quartet_count: int
    uninformative_quartet_count: int
    concordant_frequency: float
    discordant_frequency: float
    uninformative_frequency: float


@dataclass(slots=True)
class TreeSetQuartetSupportReport:
    reference_tree_path: Path
    comparison_tree_set_path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    branch_count: int
    total_quartet_count: int
    concordant_quartet_count: int
    discordant_quartet_count: int
    uninformative_quartet_count: int
    rows: list[TreeSetQuartetSupportRow]


@dataclass(slots=True)
class GeneTreeQuartetConcordanceRow:
    branch_id: str
    left_taxa: list[str]
    right_taxa: list[str]
    quartet_count_per_tree: int
    concordant_quartet_count: int
    discordant_first_quartet_count: int
    discordant_second_quartet_count: int
    uninformative_quartet_count: int
    informative_quartet_count: int
    concordance_factor: float | None
    concordant_frequency: float
    discordant_first_frequency: float
    discordant_second_frequency: float
    uninformative_frequency: float


@dataclass(slots=True)
class GeneTreeQuartetConcordanceReport:
    species_tree_path: Path
    gene_tree_set_path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    branch_count: int
    total_quartet_count: int
    concordant_quartet_count: int
    discordant_first_quartet_count: int
    discordant_second_quartet_count: int
    uninformative_quartet_count: int
    informative_quartet_count: int
    rows: list[GeneTreeQuartetConcordanceRow]


@dataclass(slots=True)
class CandidateTreeQuartetScoreReport:
    candidate_tree_path: Path
    gene_tree_set_path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    branch_count: int
    total_quartet_count: int
    concordant_quartet_count: int
    discordant_first_quartet_count: int
    discordant_second_quartet_count: int
    uninformative_quartet_count: int
    informative_quartet_count: int
    quartet_score: int
    normalized_quartet_score: float | None
    rows: list[GeneTreeQuartetConcordanceRow]


@dataclass(slots=True)
class QuartetTopologyScoreRow:
    quartet_taxa: list[str]
    first_split_taxa: list[str]
    first_split_tree_count: int
    second_split_taxa: list[str]
    second_split_tree_count: int
    third_split_taxa: list[str]
    third_split_tree_count: int
    uninformative_tree_count: int
    best_split_taxa: list[str] | None
    best_split_support_frequency: float | None
    tied_best_split_taxa: list[list[str]]


@dataclass(slots=True)
class QuartetPuzzlingAssemblyRow:
    order_index: int
    taxon_order: list[str]
    quartet_score: float
    assembled_topology_id: str
    assembled_tree_newick: str


@dataclass(slots=True)
class QuartetPuzzlingReport:
    tree_set_path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    quartet_count: int
    assembly_count: int
    unique_assembled_topology_count: int
    canonical_root_taxon: str
    canonical_rooting_strategy: str
    consensus_method: str
    consensus_threshold: float
    included_clade_count: int
    consensus_newick: str
    quartet_rows: list[QuartetTopologyScoreRow]
    assembly_rows: list[QuartetPuzzlingAssemblyRow]


@dataclass(slots=True)
class ConsensusTreeReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    consensus_method: str
    consensus_threshold: float
    included_clade_count: int
    consensus_newick: str


@dataclass(frozen=True, slots=True)
class TreeDistancePair:
    left_index: int
    right_index: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float


@dataclass(slots=True)
class TreeDistanceMatrixReport:
    path: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    shared_taxa: list[str]
    pairs: list[TreeDistancePair]
