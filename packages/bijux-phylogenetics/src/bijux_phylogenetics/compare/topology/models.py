from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.phylo.pruning import RequestedTaxaPruningReport

RobinsonFouldsMode = str
TaxonOverlapPolicy = str
BranchScoreStatus = str

_STRONG_SUPPORT_THRESHOLD = 0.9
_WEAK_SUPPORT_THRESHOLD = 0.7
_SUPPORT_DISAGREEMENT_THRESHOLD = 0.15


@dataclass(slots=True)
class RobinsonFouldsComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    rf_mode: str
    left_split_count: int
    right_split_count: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    topology_equal: bool


@dataclass(slots=True)
class TreeComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    rf_mode: str
    left_informative_clades: int
    right_informative_clades: int
    left_unrooted_splits: int
    right_unrooted_splits: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    rooted_robinson_foulds_distance: int
    rooted_normalized_robinson_foulds: float
    unrooted_robinson_foulds_distance: int
    unrooted_normalized_robinson_foulds: float
    topology_equal: bool
    same_unrooted_topology: bool
    same_taxa_different_rooting: bool
    same_topology_different_branch_lengths: bool


@dataclass(slots=True)
class CladeAgeComparisonRow:
    clade_id: str
    node_kind: str
    taxon_count: int
    descendant_taxa: list[str]
    left_age: float
    right_age: float
    age_difference: float
    absolute_age_difference: float
    unstable_age: bool


@dataclass(slots=True)
class DateAwareTreeComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    comparison_scope: str
    left_root_age: float
    right_root_age: float
    age_rmse: float
    mean_absolute_age_difference: float
    max_absolute_age_difference: float
    unstable_age_threshold: float
    unstable_clade_count: int
    matched_clade_count: int
    topology: TreeComparisonReport
    clade_rows: list[CladeAgeComparisonRow]


@dataclass(slots=True)
class DeepCoalescenceTaxonMapRow:
    gene_taxon: str
    species_taxon: str


@dataclass(slots=True)
class DeepCoalescenceBranchRow:
    species_branch: str
    branch_role: str
    descendant_species: list[str]
    lineage_count_entering: int
    coalescent_event_count: int
    lineage_count_exiting: int
    extra_lineage_count: int
    included_in_deep_coalescence_total: bool


@dataclass(slots=True)
class DeepCoalescenceReport:
    species_tree_path: Path
    gene_tree_path: Path
    observed_species_taxa: list[str]
    species_only_taxa: list[str]
    gene_tip_count: int
    deep_coalescence_total: int
    mapping_rows: list[DeepCoalescenceTaxonMapRow]
    branch_rows: list[DeepCoalescenceBranchRow]


@dataclass(slots=True)
class DuplicationLossTransferAssociationRow:
    gene_taxon: str
    species_taxon: str


@dataclass(slots=True)
class DuplicationLossTransferEventRow:
    gene_node: str
    gene_node_name: str | None
    descendant_gene_tips: list[str]
    event_type: str
    mapped_species_branch: str
    mapped_descendant_species: list[str]
    left_child_gene_node: str | None
    right_child_gene_node: str | None
    left_child_species_branch: str | None
    right_child_species_branch: str | None
    transferred_child_side: str | None
    transfer_recipient_branch: str | None
    loss_branches: list[str]
    event_cost: float


@dataclass(slots=True)
class DuplicationLossTransferReport:
    species_tree_path: Path
    gene_tree_path: Path
    taxon_map_path: Path | None
    duplication_cost: float
    loss_cost: float
    transfer_cost: float
    observed_species_taxa: list[str]
    species_only_taxa: list[str]
    gene_tip_count: int
    root_mapping_branch: str
    reconciliation_score: float
    duplication_event_count: int
    loss_event_count: int
    transfer_event_count: int
    speciation_event_count: int
    mapping_rows: list[DuplicationLossTransferAssociationRow]
    event_rows: list[DuplicationLossTransferEventRow]


@dataclass(slots=True)
class SharedTaxaPruningReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    left_pruning: RequestedTaxaPruningReport
    right_pruning: RequestedTaxaPruningReport
    post_pruning_comparison: TreeComparisonReport


@dataclass(slots=True)
class AgreementSubtreeCandidateRow:
    candidate_index: int
    retained_taxon_count: int
    retained_taxa: list[str]
    removed_taxa: list[str]
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    topology_equal: bool


@dataclass(slots=True)
class AgreementSubtreePruningReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    rf_mode: str
    search_strategy: str
    possible_retained_subset_count: int
    evaluated_candidate_count: int
    retained_taxa: list[str]
    agreement_removed_taxa: list[str]
    left_pruning: RequestedTaxaPruningReport
    right_pruning: RequestedTaxaPruningReport
    post_pruning_comparison: TreeComparisonReport
    candidate_rows: list[AgreementSubtreeCandidateRow]


@dataclass(slots=True)
class MaximumAgreementSubtreeSearchRow:
    evaluation_index: int
    step_index: int
    retained_taxon_count: int
    retained_taxa: list[str]
    removed_taxa: list[str]
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    topology_equal: bool
    selected_for_next_step: bool


@dataclass(slots=True)
class MaximumAgreementSubtreeApproximationReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    rf_mode: str
    search_strategy: str
    selection_objective: str
    approximation_status: str
    possible_retained_subset_count: int
    max_evaluated_candidate_count: int
    evaluated_candidate_count: int
    retained_taxa: list[str]
    approximation_removed_taxa: list[str]
    left_pruning: RequestedTaxaPruningReport
    right_pruning: RequestedTaxaPruningReport
    post_pruning_comparison: TreeComparisonReport
    search_rows: list[MaximumAgreementSubtreeSearchRow]


@dataclass(slots=True)
class CladeSupportPair:
    split_id: str
    left_support: float | None
    right_support: float | None
    left_support_fraction: float | None
    right_support_fraction: float | None
    support_fraction_delta: float | None
    support_disagreement: bool


@dataclass(slots=True)
class SupportConflictRow:
    split_id: str
    comparison_status: str
    left_present: bool
    right_present: bool
    left_support: float | None
    right_support: float | None
    left_support_fraction: float | None
    right_support_fraction: float | None
    strongest_support_fraction: float | None
    support_strength: str
    conflict_classification: str
    detail: str


@dataclass(slots=True)
class SupportComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    strong_support_threshold: float
    weak_support_threshold: float
    support_disagreement_threshold: float
    shared_clades: list[CladeSupportPair]
    conflicting_clades: list[SupportConflictRow]


@dataclass(slots=True)
class BranchLengthPair:
    split_id: str
    left_length: float | None
    right_length: float | None
    delta: float | None
    ratio: float | None


@dataclass(slots=True)
class BranchScoreSplit:
    split_id: str
    comparison_status: str
    left_length: float | None
    right_length: float | None
    branch_score_difference: float | None
    squared_difference: float | None


@dataclass(slots=True)
class BranchScoreComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    same_taxon_set: bool
    branch_score_distance: float | None
    split_count: int
    shared_split_count: int
    left_only_split_count: int
    right_only_split_count: int
    missing_length_split_count: int
    splits: list[BranchScoreSplit]


@dataclass(slots=True)
class BranchLengthComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    same_taxon_set: bool
    shared_splits: list[BranchLengthPair]
    branch_score: BranchScoreComparisonReport


@dataclass(slots=True)
class CladeSetComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    shared_clades: list[str]
    left_only_clades: list[str]
    right_only_clades: list[str]


@dataclass(slots=True)
class CladeOverlapObservation:
    tree_path: Path
    present: bool
    support: float | None


@dataclass(slots=True)
class CladeOverlapRow:
    clade_id: str
    present_in_all_trees: bool
    present_tree_count: int
    absent_tree_count: int
    observations: list[CladeOverlapObservation]


@dataclass(slots=True)
class TreeCladeOverlapSummary:
    tree_path: Path
    clade_count: int
    support_clade_count: int
    unique_clades: list[str]
    excluded_taxa: list[str]


@dataclass(slots=True)
class CladeOverlapComparisonReport:
    tree_paths: list[Path]
    shared_taxa: list[str]
    shared_clades: list[str]
    conflicting_clades: list[str]
    tree_summaries: list[TreeCladeOverlapSummary]
    clade_rows: list[CladeOverlapRow]


@dataclass(slots=True)
class CladeChangeReport:
    left_path: Path
    right_path: Path
    lost_clades: list[str]
    gained_clades: list[str]


@dataclass(slots=True)
class InMemoryTopologyComparison:
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    rf_mode: str
    left_informative_clades: int
    right_informative_clades: int
    left_unrooted_splits: int
    right_unrooted_splits: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    rooted_robinson_foulds_distance: int
    rooted_normalized_robinson_foulds: float
    unrooted_robinson_foulds_distance: int
    unrooted_normalized_robinson_foulds: float
    topology_equal: bool
    same_unrooted_topology: bool
    same_taxa_different_rooting: bool


@dataclass(slots=True)
class InMemoryBranchLengthComparison:
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    taxon_overlap_policy: str
    shared_splits: list[BranchLengthPair]
    branch_score_distance: float | None
    branch_score_splits: list[BranchScoreSplit]
    missing_length_split_count: int
