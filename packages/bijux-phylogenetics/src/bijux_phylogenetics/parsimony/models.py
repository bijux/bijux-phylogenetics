from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class FitchCharacterMatrix:
    """Taxon-by-character matrix normalized for unordered Fitch scoring."""

    matrix_path: Path | None
    taxon_column: str
    character_ids: list[str]
    states_by_taxon: dict[str, dict[str, str]]

    @property
    def taxon_count(self) -> int:
        return len(self.states_by_taxon)

    @property
    def character_count(self) -> int:
        return len(self.character_ids)


@dataclass(frozen=True, slots=True)
class FitchCharacterScore:
    """Per-character unordered Fitch tree-length row."""

    character_id: str
    step_count: int
    observed_states: list[str]
    character_weight: float
    weighted_score: float


@dataclass(frozen=True, slots=True)
class FitchNodeStateSet:
    """Per-node unordered Fitch candidate-state row for one character."""

    character_id: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    state_set: list[str]
    is_tip: bool
    observed_state: str | None


@dataclass(frozen=True, slots=True)
class FitchScoreReport:
    """Complete unordered Fitch scoring report over one tree and matrix."""

    algorithm: str
    tree_path: Path | None
    matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    total_steps: int
    weights_path: Path | None
    total_weighted_score: float
    step_rows: list[FitchCharacterScore]
    node_state_rows: list[FitchNodeStateSet]


@dataclass(frozen=True, slots=True)
class ParsimonyPlacementAlternativeRow:
    """One query placement scored on one candidate edge by additional Fitch steps."""

    query_id: str
    placement_rank: int
    edge_id: str
    child_name: str | None
    descendant_taxa: list[str]
    total_steps: int
    additional_steps: int
    total_weighted_score: float
    additional_weighted_score: float
    is_equally_best: bool
    placed_tree_newick: str


@dataclass(frozen=True, slots=True)
class ParsimonyPlacementQuerySummary:
    """Best-placement summary for one parsimony query taxon."""

    query_id: str
    character_count: int
    best_edge_id: str
    best_child_name: str | None
    best_descendant_taxa: list[str]
    best_total_steps: int
    best_additional_steps: int
    best_total_weighted_score: float
    best_additional_weighted_score: float
    candidate_placement_count: int
    equally_best_placement_count: int
    selected_best_tree_newick: str


@dataclass(frozen=True, slots=True)
class ParsimonyPlacementReport:
    """Exact unordered-Fitch placement report for one or more query taxa."""

    algorithm: str
    method: str
    tree_path: Path | None
    matrix_path: Path | None
    query_matrix_path: Path | None
    taxon_column: str
    reference_taxon_count: int
    character_count: int
    edge_count: int
    query_count: int
    reference_total_steps: int
    weights_path: Path | None
    reference_total_weighted_score: float
    query_summaries: list[ParsimonyPlacementQuerySummary]
    alternative_rows: list[ParsimonyPlacementAlternativeRow]


ParsimonyCharacterMatrix = FitchCharacterMatrix


@dataclass(frozen=True, slots=True)
class ParsimonyCharacterWeights:
    """Validated per-character weights for one governed parsimony matrix."""

    weights_path: Path | None
    weights_by_character: dict[str, float]


@dataclass(frozen=True, slots=True)
class WagnerCharacterScore:
    """Per-character ordered Wagner weighted tree-length row."""

    character_id: str
    weighted_step_count: int
    observed_states: list[str]
    state_order: list[str]
    optimal_root_states: list[str]
    character_weight: float
    weighted_score: float


@dataclass(frozen=True, slots=True)
class WagnerNodeCost:
    """Per-node ordered Wagner cost row for one candidate state."""

    character_id: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    state: str
    cost: int
    is_optimal_state: bool


@dataclass(frozen=True, slots=True)
class WagnerScoreReport:
    """Complete ordered Wagner scoring report over one tree and matrix."""

    algorithm: str
    tree_path: Path | None
    matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    total_cost: int
    weights_path: Path | None
    total_weighted_score: float
    step_rows: list[WagnerCharacterScore]
    node_cost_rows: list[WagnerNodeCost]


@dataclass(frozen=True, slots=True)
class SankoffCostMatrix:
    """Validated state-to-state Sankoff transition-cost matrix."""

    matrix_path: Path | None
    states: list[str]
    costs: dict[str, dict[str, float]]
    validation_warnings: list[SankoffCostMatrixWarning]


@dataclass(frozen=True, slots=True)
class SankoffCostMatrixWarning:
    """Structured nonfatal warning surfaced while validating one Sankoff cost matrix."""

    code: str
    message: str
    details: dict[str, object]


@dataclass(frozen=True, slots=True)
class SankoffCharacterScore:
    """Per-character Sankoff minimum-cost row."""

    character_id: str
    minimum_cost: float
    observed_states: list[str]
    matrix_states: list[str]
    character_weight: float
    weighted_score: float


@dataclass(frozen=True, slots=True)
class SankoffNodeCost:
    """Per-node Sankoff cost-vector row for one candidate state."""

    character_id: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    state: str
    cost: float
    is_optimal_state: bool


@dataclass(frozen=True, slots=True)
class SankoffNodeSelection:
    """Per-node Sankoff optimal-state selection row."""

    character_id: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    optimal_states: list[str]
    tie_states: list[str]


@dataclass(frozen=True, slots=True)
class SankoffScoreReport:
    """Complete Sankoff scoring report over one tree, matrix, and cost matrix."""

    algorithm: str
    tree_path: Path | None
    matrix_path: Path | None
    cost_matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    total_cost: float
    weights_path: Path | None
    total_weighted_score: float
    validation_warnings: list[SankoffCostMatrixWarning]
    step_rows: list[SankoffCharacterScore]
    node_cost_rows: list[SankoffNodeCost]
    node_selection_rows: list[SankoffNodeSelection]


@dataclass(frozen=True, slots=True)
class DolloCharacterScore:
    """Per-character Dollo summary row."""

    character_id: str
    step_count: int
    derived_taxon_count: int
    gain_node: str | None
    gain_node_name: str | None
    gain_descendant_taxa: list[str]
    total_losses: int
    impossible_state_warning: str | None
    character_weight: float
    weighted_score: float


@dataclass(frozen=True, slots=True)
class DolloBranchChange:
    """Per-branch Dollo change row."""

    character_id: str
    change_kind: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]


@dataclass(frozen=True, slots=True)
class DolloScoreReport:
    """Complete Dollo scoring report over one tree and binary matrix."""

    algorithm: str
    tree_path: Path | None
    matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    total_gains: int
    total_losses: int
    weights_path: Path | None
    total_weighted_score: float
    step_rows: list[DolloCharacterScore]
    branch_change_rows: list[DolloBranchChange]


@dataclass(frozen=True, slots=True)
class CaminSokalCharacterScore:
    """Per-character irreversible Camin-Sokal summary row."""

    character_id: str
    derived_taxon_count: int
    gain_count: int
    root_state: str
    character_weight: float
    weighted_score: float


@dataclass(frozen=True, slots=True)
class CaminSokalBranchChange:
    """Per-branch irreversible Camin-Sokal change row."""

    character_id: str
    change_kind: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]


@dataclass(frozen=True, slots=True)
class CaminSokalScoreReport:
    """Complete irreversible Camin-Sokal scoring report over one tree and binary matrix."""

    algorithm: str
    tree_path: Path | None
    matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    root_state: str
    total_gains: int
    weights_path: Path | None
    total_weighted_score: float
    step_rows: list[CaminSokalCharacterScore]
    branch_change_rows: list[CaminSokalBranchChange]


@dataclass(frozen=True, slots=True)
class ParsimonyReconstructionCharacterScore:
    """Per-character unordered parsimony reconstruction summary row."""

    character_id: str
    step_count: int
    observed_states: list[str]
    root_state: str
    character_weight: float
    weighted_score: float


@dataclass(frozen=True, slots=True)
class ParsimonyReconstructionNodeState:
    """Per-node resolved parsimony state row for one character."""

    character_id: str
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    resolved_state: str
    is_tip: bool
    observed_state: str | None


@dataclass(frozen=True, slots=True)
class ParsimonyAncestralState:
    """Per-node ancestral-state export row for one parsimony reconstruction."""

    node_id: str
    clade_id: str
    character_id: str
    possible_states: list[str]
    chosen_state: str
    method: str
    ambiguous: bool


@dataclass(frozen=True, slots=True)
class ParsimonyReconstructionBranchChange:
    """Per-branch resolved parsimony change row for one character."""

    character_id: str
    branch_id: str
    parent_node: str
    parent_state: str
    child_node: str
    child_node_name: str | None
    child_descendant_taxa: list[str]
    change_from: str
    change_to: str
    ambiguous: bool

    @property
    def node(self) -> str:
        """Backward-compatible alias for the child node branch endpoint."""
        return self.child_node

    @property
    def node_name(self) -> str | None:
        """Backward-compatible alias for the child node label."""
        return self.child_node_name

    @property
    def descendant_taxa(self) -> list[str]:
        """Backward-compatible alias for the child node descendant taxa."""
        return self.child_descendant_taxa


@dataclass(frozen=True, slots=True)
class ParsimonyReconstructionReport:
    """Complete unordered parsimony reconstruction report over one tree and matrix."""

    algorithm: str
    tree_path: Path | None
    matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    total_steps: int
    weights_path: Path | None
    total_weighted_score: float
    step_rows: list[ParsimonyReconstructionCharacterScore]
    node_state_rows: list[ParsimonyReconstructionNodeState]
    ancestral_state_rows: list[ParsimonyAncestralState]
    branch_change_rows: list[ParsimonyReconstructionBranchChange]


@dataclass(frozen=True, slots=True)
class ParsimonyTreeLengthCharacterScore:
    """Per-character raw and weighted tree-length row."""

    character_id: str
    raw_score: float
    character_weight: float
    weighted_score: float


@dataclass(frozen=True, slots=True)
class ParsimonyTreeLengthReport:
    """Complete tree-length summary over one matrix scored by one parsimony method."""

    algorithm: str
    method: str
    tree_path: Path | None
    matrix_path: Path | None
    cost_matrix_path: Path | None
    weights_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    raw_total_score: float
    total_score: float
    step_rows: list[ParsimonyTreeLengthCharacterScore]


@dataclass(frozen=True, slots=True)
class ParsimonyConsistencyCharacterIndex:
    """Per-character consistency-index row."""

    character_id: str
    character_kind: str
    observed_states: list[str]
    minimum_possible_steps: float
    observed_steps: float
    consistency_index: float | None
    undefined_reason: str | None


@dataclass(frozen=True, slots=True)
class ParsimonyConsistencyIndexReport:
    """Complete consistency-index summary over one matrix and parsimony method."""

    algorithm: str
    method: str
    tree_path: Path | None
    matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    included_character_count: int
    excluded_character_count: int
    minimum_possible_steps_total: float
    observed_steps_total: float
    consistency_index: float | None
    undefined_reason: str | None
    character_rows: list[ParsimonyConsistencyCharacterIndex]


@dataclass(frozen=True, slots=True)
class ParsimonyRetentionCharacterIndex:
    """Per-character retention-index row."""

    character_id: str
    character_kind: str
    observed_states: list[str]
    minimum_possible_steps: float
    maximum_possible_steps: float
    observed_steps: float
    retention_index: float | None
    undefined_reason: str | None


@dataclass(frozen=True, slots=True)
class ParsimonyRetentionIndexReport:
    """Complete retention-index summary over one matrix and parsimony method."""

    algorithm: str
    method: str
    tree_path: Path | None
    matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    included_character_count: int
    excluded_character_count: int
    minimum_possible_steps_total: float
    maximum_possible_steps_total: float
    observed_steps_total: float
    retention_index: float | None
    undefined_reason: str | None
    character_rows: list[ParsimonyRetentionCharacterIndex]


@dataclass(frozen=True, slots=True)
class ParsimonyRescaledConsistencyCharacterIndex:
    """Per-character rescaled consistency-index row."""

    character_id: str
    ci: float | None
    ri: float | None
    rc: float | None
    undefined_reason: str | None


@dataclass(frozen=True, slots=True)
class ParsimonyRescaledConsistencyIndexReport:
    """Complete rescaled consistency-index summary over one matrix and method."""

    algorithm: str
    method: str
    tree_path: Path | None
    matrix_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    ci: float | None
    ri: float | None
    rc: float | None
    undefined_reason: str | None
    character_rows: list[ParsimonyRescaledConsistencyCharacterIndex]


@dataclass(frozen=True, slots=True)
class ParsimonyBootstrapReplicate:
    """One exact-search bootstrap replicate and its sampled characters."""

    replicate_index: int
    sampled_character_ids: list[str]
    best_score: float
    optimal_tree_count: int
    tree_newick: str


@dataclass(frozen=True, slots=True)
class ParsimonyBootstrapCladeSupport:
    """One clade-identity bootstrap support row on the reference tree."""

    branch_id: str
    node_name: str | None
    descendant_taxa: list[str]
    supporting_tree_count: int
    clade_frequency: float
    support_percent: float


@dataclass(frozen=True, slots=True)
class ParsimonyBootstrapReport:
    """Complete exact-search parsimony bootstrap report over one character matrix."""

    algorithm: str
    method: str
    matrix_path: Path | None
    cost_matrix_path: Path | None
    weights_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    replicate_count: int
    random_seed: int
    candidate_tree_count: int
    max_exact_taxa: int
    reference_score: float
    reference_optimal_tree_count: int
    reference_tree_newick: str
    replicate_rows: list[ParsimonyBootstrapReplicate]
    clade_support_rows: list[ParsimonyBootstrapCladeSupport]


@dataclass(frozen=True, slots=True)
class ParsimonyJackknifeReplicate:
    """One exact-search jackknife replicate and its retained characters."""

    replicate_index: int
    retained_character_count: int
    retained_character_ids: list[str]
    best_score: float
    optimal_tree_count: int
    tree_newick: str


@dataclass(frozen=True, slots=True)
class ParsimonyJackknifeCladeSupport:
    """One clade-identity jackknife support row on the reference tree."""

    branch_id: str
    node_name: str | None
    descendant_taxa: list[str]
    supporting_tree_count: int
    clade_frequency: float
    support_percent: float


@dataclass(frozen=True, slots=True)
class ParsimonyJackknifeReport:
    """Complete exact-search parsimony jackknife report over one character matrix."""

    algorithm: str
    method: str
    matrix_path: Path | None
    cost_matrix_path: Path | None
    weights_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    replicate_count: int
    random_seed: int
    retain_probability: float
    candidate_tree_count: int
    max_exact_taxa: int
    reference_score: float
    reference_optimal_tree_count: int
    reference_tree_newick: str
    replicate_rows: list[ParsimonyJackknifeReplicate]
    clade_support_rows: list[ParsimonyJackknifeCladeSupport]


@dataclass(frozen=True, slots=True)
class ParsimonyEqualBestTree:
    """One retained equally optimal tree from a deterministic exact parsimony search."""

    tree_index: int
    total_score: float
    tree_newick: str


@dataclass(frozen=True, slots=True)
class ParsimonyConsensusSummary:
    """One consensus summary over a retained equal-best tree set."""

    consensus_method: str
    consensus_threshold: float
    tree_count: int
    included_clade_count: int
    consensus_newick: str


@dataclass(frozen=True, slots=True)
class ParsimonyEqualBestConsensusReport:
    """Exact parsimony search summary over the retained equally optimal tree set."""

    algorithm: str
    method: str
    matrix_path: Path | None
    cost_matrix_path: Path | None
    weights_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    candidate_tree_count: int
    max_exact_taxa: int
    max_retained_equal_best_trees: int
    best_score: float
    equal_best_tree_count: int
    retained_equal_best_tree_count: int
    retained_all_equal_best_trees: bool
    strict_consensus: ParsimonyConsensusSummary | None
    majority_consensus: ParsimonyConsensusSummary | None
    equal_best_tree_rows: list[ParsimonyEqualBestTree]


@dataclass(frozen=True, slots=True)
class ParsimonyBremerSupportRow:
    """One rooted reference-tree clade with its exact Bremer decay summary."""

    branch_id: str
    node_name: str | None
    descendant_taxa: list[str]
    shortest_lacking_score: float
    decay_index: float
    shortest_lacking_tree_count: int
    shortest_lacking_tree_newick: str


@dataclass(frozen=True, slots=True)
class ParsimonyBremerSupportReport:
    """Exact small-taxon Bremer support summary over one reference tree and matrix."""

    algorithm: str
    method: str
    tree_path: Path | None
    matrix_path: Path | None
    cost_matrix_path: Path | None
    weights_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    candidate_tree_count: int
    max_exact_taxa: int
    reference_tree_newick: str
    reference_tree_score: float
    optimal_score: float
    optimal_tree_count: int
    optimal_tree_newick: str
    reference_tree_score_delta_from_optimal: float
    reference_tree_is_optimal: bool
    bremer_rows: list[ParsimonyBremerSupportRow]


@dataclass(frozen=True, slots=True)
class ParsimonyNniTraceRow:
    """One deterministic event row in a rooted parsimony NNI search trace."""

    event_index: int
    event_kind: str
    iteration: int
    score_before: float | None
    score_after: float
    score_delta: float | None
    tree_before_newick: str | None
    tree_after_newick: str
    pivot_branch_id: str | None
    sibling_clade_id: str | None
    exchanged_clade_id: str | None
    stopping_reason: str | None


@dataclass(frozen=True, slots=True)
class ParsimonyNniSearchReport:
    """Complete rooted NNI parsimony hill-climb report over one starting tree."""

    algorithm: str
    method: str
    tree_path: Path | None
    matrix_path: Path | None
    cost_matrix_path: Path | None
    weights_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    start_tree_newick: str
    start_score: float
    final_tree_newick: str
    final_score: float
    accepted_move_count: int
    evaluated_neighbor_count: int
    stopping_reason: str
    trace_rows: list[ParsimonyNniTraceRow]


@dataclass(frozen=True, slots=True)
class ParsimonySprTraceRow:
    """One deterministic event row in a rooted parsimony SPR search trace."""

    event_index: int
    event_kind: str
    iteration: int
    score_before: float | None
    score_after: float
    score_delta: float | None
    tree_before_newick: str | None
    tree_after_newick: str
    pruned_clade_id: str | None
    regraft_target_branch_id: str | None
    stopping_reason: str | None


@dataclass(frozen=True, slots=True)
class ParsimonySprSearchReport:
    """Complete rooted SPR parsimony hill-climb report over one starting tree."""

    algorithm: str
    method: str
    tree_path: Path | None
    matrix_path: Path | None
    cost_matrix_path: Path | None
    weights_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    start_tree_newick: str
    start_score: float
    final_tree_newick: str
    final_score: float
    accepted_move_count: int
    evaluated_neighbor_count: int
    stopping_reason: str
    trace_rows: list[ParsimonySprTraceRow]


@dataclass(frozen=True, slots=True)
class ParsimonyRatchetCycle:
    """One deterministic ratchet cycle over perturbed then restored character weights."""

    cycle_index: int
    start_score: float
    start_tree_newick: str
    perturbed_character_ids: list[str]
    perturbation_factor: float
    perturbed_score: float
    perturbed_tree_newick: str
    perturbed_accepted_move_count: int
    restored_score: float
    restored_tree_newick: str
    restored_accepted_move_count: int
    best_score_after_cycle: float
    best_tree_after_cycle: str
    best_tree_improved: bool


@dataclass(frozen=True, slots=True)
class ParsimonyRatchetBestTreeHistory:
    """One deterministic best-tree history row across ratchet cycles."""

    history_index: int
    cycle_index: int
    best_score: float
    best_tree_newick: str


@dataclass(frozen=True, slots=True)
class ParsimonyRatchetReport:
    """Complete parsimony ratchet report over one starting tree and character matrix."""

    algorithm: str
    method: str
    tree_path: Path | None
    matrix_path: Path | None
    cost_matrix_path: Path | None
    weights_path: Path | None
    taxon_column: str
    taxon_count: int
    character_count: int
    cycle_count: int
    random_seed: int
    perturbed_character_count: int
    perturbation_factor: float
    start_tree_newick: str
    start_score: float
    final_tree_newick: str
    final_score: float
    best_tree_newick: str
    best_score: float
    cycle_rows: list[ParsimonyRatchetCycle]
    best_tree_history_rows: list[ParsimonyRatchetBestTreeHistory]
