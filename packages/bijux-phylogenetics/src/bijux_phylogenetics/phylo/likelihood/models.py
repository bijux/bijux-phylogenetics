from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Jc69BranchLengthOptimizationStep:
    """One bounded branch-length update inside a fixed-topology JC69 search."""

    optimization_pass: int
    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    starting_branch_length: float
    optimized_branch_length: float
    starting_log_likelihood: float
    optimized_log_likelihood: float
    accepted: bool


@dataclass(slots=True)
class Jc69TreeLikelihoodReport:
    """Native JC69 likelihood report for one fixed topology and alignment."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    state_count: int
    observation_policy: str
    log_likelihood: float


@dataclass(slots=True)
class Jc69BranchLengthOptimizationReport:
    """Fixed-topology JC69 branch-length optimization summary."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    branch_count: int
    initial_tree_newick: str
    optimized_tree_newick: str
    initial_log_likelihood: float
    optimized_log_likelihood: float
    optimization_pass_count: int
    function_evaluation_count: int
    converged: bool
    lower_branch_length_bound: float
    upper_branch_length_bound: float
    steps: list[Jc69BranchLengthOptimizationStep]


@dataclass(frozen=True, slots=True)
class StrictClockBranchRow:
    """One branch duration and its optimized strict-clock substitution length."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    time_duration: float
    optimized_branch_length: float
    optimized_clock_rate: float


@dataclass(slots=True)
class StrictClockLikelihoodReport:
    """One fixed-topology strict-clock likelihood fit on a time-scaled tree."""

    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    branch_count: int
    compression_used: bool
    time_tree_newick: str
    scaled_tree_newick: str
    initial_clock_rate: float
    optimized_clock_rate: float
    initial_log_likelihood: float
    optimized_log_likelihood: float
    parameter_count: int
    aic: float
    function_evaluation_count: int
    converged: bool
    lower_clock_rate_bound: float
    upper_clock_rate_bound: float
    branch_rows: list[StrictClockBranchRow]
    site_log_likelihoods: list[SiteLogLikelihoodRow]
    branch_likelihood_diagnostics: list[BranchLikelihoodDiagnosticRow]


@dataclass(frozen=True, slots=True)
class LocalClockRegimeRow:
    """One fitted local-clock regime and the branches it owns."""

    regime_id: str
    target_kind: str
    target_label: str | None
    descendant_taxa: list[str]
    node_id: str | None
    branch_count: int
    optimized_clock_rate: float


@dataclass(frozen=True, slots=True)
class LocalClockBranchRow:
    """One branch assignment under a fitted local-clock likelihood model."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    regime_id: str
    target_kind: str
    time_duration: float
    optimized_branch_length: float
    optimized_clock_rate: float


@dataclass(slots=True)
class LocalClockLikelihoodReport:
    """One fixed-topology local-clock JC69 likelihood fit on a time-scaled tree."""

    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    branch_count: int
    regime_count: int
    compression_used: bool
    time_tree_newick: str
    scaled_tree_newick: str
    strict_clock_rate: float
    strict_clock_log_likelihood: float
    strict_clock_aic: float
    initial_clock_rate: float
    optimized_log_likelihood: float
    parameter_count: int
    aic: float
    aic_delta_vs_strict_clock: float
    preferred_model_by_aic: str
    function_evaluation_count: int
    optimization_pass_count: int
    converged: bool
    lower_clock_rate_bound: float
    upper_clock_rate_bound: float
    branch_rows: list[LocalClockBranchRow]
    regime_rows: list[LocalClockRegimeRow]
    site_log_likelihoods: list[SiteLogLikelihoodRow]
    branch_likelihood_diagnostics: list[BranchLikelihoodDiagnosticRow]


@dataclass(slots=True)
class K80TreeLikelihoodReport:
    """Native K80 likelihood report for one fixed topology and alignment."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    state_count: int
    observation_policy: str
    kappa: float
    log_likelihood: float


@dataclass(slots=True)
class K80KappaOptimizationReport:
    """Fixed-topology K80 kappa optimization summary."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    tree_newick: str
    initial_kappa: float
    optimized_kappa: float
    initial_log_likelihood: float
    optimized_log_likelihood: float
    function_evaluation_count: int
    converged: bool
    lower_kappa_bound: float
    upper_kappa_bound: float


@dataclass(slots=True)
class F81TreeLikelihoodReport:
    """Native F81 likelihood report for one fixed topology and alignment."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    state_count: int
    observation_policy: str
    base_frequency_source: str
    base_frequency_a: float
    base_frequency_c: float
    base_frequency_g: float
    base_frequency_t: float
    parameter_count: int
    log_likelihood: float
    aic: float


@dataclass(slots=True)
class Hky85TreeLikelihoodReport:
    """Native HKY85 likelihood report for one fixed topology and alignment."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    state_count: int
    observation_policy: str
    base_frequency_source: str
    base_frequency_a: float
    base_frequency_c: float
    base_frequency_g: float
    base_frequency_t: float
    kappa: float
    parameter_count: int
    log_likelihood: float
    aic: float


@dataclass(slots=True)
class Hky85KappaOptimizationReport:
    """Fixed-topology HKY85 kappa optimization summary."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    tree_newick: str
    base_frequency_source: str
    base_frequency_a: float
    base_frequency_c: float
    base_frequency_g: float
    base_frequency_t: float
    initial_kappa: float
    optimized_kappa: float
    parameter_count: int
    initial_log_likelihood: float
    optimized_log_likelihood: float
    initial_aic: float
    optimized_aic: float
    function_evaluation_count: int
    converged: bool
    lower_kappa_bound: float
    upper_kappa_bound: float


@dataclass(slots=True)
class GtrTreeLikelihoodReport:
    """Native GTR likelihood report for one fixed topology and alignment."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    state_count: int
    observation_policy: str
    base_frequency_source: str
    base_frequency_a: float
    base_frequency_c: float
    base_frequency_g: float
    base_frequency_t: float
    exchangeability_anchor: str
    exchangeability_ac: float
    exchangeability_ag: float
    exchangeability_at: float
    exchangeability_cg: float
    exchangeability_ct: float
    exchangeability_gt: float
    parameter_count: int
    log_likelihood: float
    aic: float


@dataclass(slots=True)
class GtrExchangeabilityOptimizationReport:
    """Fixed-topology GTR exchangeability optimization summary."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    tree_newick: str
    base_frequency_source: str
    base_frequency_a: float
    base_frequency_c: float
    base_frequency_g: float
    base_frequency_t: float
    exchangeability_anchor: str
    exchangeability_ac: float
    exchangeability_ag: float
    exchangeability_at: float
    exchangeability_cg: float
    exchangeability_ct: float
    exchangeability_gt: float
    parameter_count: int
    initial_log_likelihood: float
    optimized_log_likelihood: float
    initial_aic: float
    optimized_aic: float
    function_evaluation_count: int
    optimization_pass_count: int
    converged: bool
    lower_exchangeability_bound: float
    upper_exchangeability_bound: float


@dataclass(slots=True)
class SubstitutionParameterOptimizationRow:
    """One optimized substitution parameter with explicit search bounds."""

    parameter_name: str
    initial_value: float
    optimized_value: float
    lower_bound: float
    upper_bound: float
    hit_lower_bound: bool
    hit_upper_bound: bool


@dataclass(slots=True)
class LikelihoodOptimizationBoundaryWarning:
    """One explicit boundary warning emitted by a likelihood optimization surface."""

    warning_kind: str
    affected_parameter: str | None
    boundary_side: str
    observed_value: float
    lower_bound: float | None
    upper_bound: float | None
    affected_branch_id: str | None
    affected_branch_clade_id: str | None
    message: str


@dataclass(slots=True)
class NucleotideSubstitutionParameterOptimizationReport:
    """Fixed-topology nucleotide substitution-parameter optimization summary."""

    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    tree_newick: str
    parameter_count: int
    base_frequency_source: str | None
    base_frequency_a: float | None
    base_frequency_c: float | None
    base_frequency_g: float | None
    base_frequency_t: float | None
    fixed_parameter_values: dict[str, float]
    parameter_rows: list[SubstitutionParameterOptimizationRow]
    initial_log_likelihood: float
    optimized_log_likelihood: float
    initial_aic: float
    optimized_aic: float
    function_evaluation_count: int
    optimization_pass_count: int
    converged: bool
    boundary_warnings: list[LikelihoodOptimizationBoundaryWarning]
    warnings: list[str]


@dataclass(slots=True)
class NestedLikelihoodRatioModelFit:
    """One fitted nucleotide model summary used in one nested likelihood-ratio test."""

    model_name: str
    parameter_count: int
    log_likelihood: float
    aic: float
    parameter_values: dict[str, float]
    warnings: list[str]


@dataclass(slots=True)
class NestedLikelihoodRatioTestReport:
    """One declared nested fixed-topology nucleotide likelihood-ratio test report."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    tree_newick: str
    null_fit: NestedLikelihoodRatioModelFit
    alternative_fit: NestedLikelihoodRatioModelFit
    likelihood_ratio_statistic: float
    degrees_of_freedom: int
    p_value: float
    p_value_method: str
    boundary_caveat: str
    warnings: list[str]


@dataclass(slots=True)
class SubstitutionModelSelectionRow:
    """One fixed-topology substitution-model selection row for one candidate surface."""

    model_name: str
    base_model_name: str
    rate_heterogeneity_model: str
    fit_succeeded: bool
    parameter_count: int | None
    log_likelihood: float | None
    aic: float | None
    aicc: float | None
    bic: float | None
    delta_aic: float | None
    akaike_weight: float | None
    rank: int | None
    comparable_on_aic: bool
    comparable_on_aicc: bool
    comparable_on_bic: bool
    selected_by_aic: bool
    selected_by_aicc: bool
    selected_by_bic: bool
    parameter_values: dict[str, float]
    warnings: list[str]


@dataclass(slots=True)
class SubstitutionModelSelectionReport:
    """Ranked fixed-topology substitution-model selection table for one alignment."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    tree_newick: str
    rows: list[SubstitutionModelSelectionRow]
    best_model_aic: str | None
    best_model_aicc: str | None
    best_model_bic: str | None
    warnings: list[str]


@dataclass(slots=True)
class NucleotideLikelihoodNniTraceRow:
    """One deterministic event row in a rooted nucleotide likelihood NNI search trace."""

    event_index: int
    event_kind: str
    iteration: int
    move_type: str
    candidate_topology_fingerprint: str
    log_likelihood_before: float | None
    log_likelihood_after: float
    log_likelihood_delta: float | None
    accepted_move: bool
    trace_reason: str
    tree_before_newick: str | None
    tree_after_newick: str
    pivot_branch_id: str | None
    sibling_clade_id: str | None
    exchanged_clade_id: str | None
    branch_reoptimization_policy: str
    branch_reoptimization_scope: str
    optimized_branch_count: int
    optimized_branch_clade_ids: list[str]
    branch_reoptimization_converged: bool | None
    branch_optimization_pass_count: int
    branch_function_evaluation_count: int
    boundary_warning_messages: list[str]
    stopping_reason: str | None


@dataclass(slots=True)
class NucleotideLikelihoodSearchConvergenceDecision:
    """Resolved stop-or-continue decision for one native likelihood tree search iteration."""

    should_stop: bool
    stopping_reason: str | None


@dataclass(slots=True)
class NucleotideLikelihoodEqualBestTreeRow:
    """One retained equal-best search tree within one native likelihood search run."""

    retained_rank: int
    topology_fingerprint: str
    tree_newick: str
    log_likelihood: float


@dataclass(slots=True)
class NucleotideLikelihoodEqualBestTreeReport:
    """Retained equal-best search trees and strict consensus for one search run."""

    likelihood_tolerance: float
    retention_cap: int
    retained_tree_count: int
    omitted_tree_count: int
    best_log_likelihood: float
    consensus_method: str
    consensus_newick: str
    rows: list[NucleotideLikelihoodEqualBestTreeRow]


@dataclass(slots=True)
class NucleotideLikelihoodNniCandidateRow:
    """One evaluated rooted NNI neighbor inside one best-improvement iteration."""

    iteration: int
    candidate_order: int
    pivot_branch_id: str
    sibling_clade_id: str
    exchanged_clade_id: str
    candidate_tree_newick: str
    log_likelihood: float
    log_likelihood_delta: float
    improving_move: bool
    selected_best_move: bool
    branch_reoptimization_scope: str
    optimized_branch_count: int
    optimized_branch_clade_ids: list[str]
    branch_reoptimization_converged: bool
    branch_optimization_pass_count: int
    branch_function_evaluation_count: int
    boundary_warning_messages: list[str]


@dataclass(slots=True)
class NucleotideLikelihoodNniSearchReport:
    """Complete rooted nucleotide likelihood NNI hill-climb report."""

    algorithm: str
    model_name: str
    improvement_policy: str
    tree_path: str | None
    alignment_path: str | None
    taxon_count: int
    site_count: int
    pattern_count: int
    input_tree_newick: str
    start_tree_newick: str
    start_log_likelihood: float
    final_tree_newick: str
    final_log_likelihood: float
    accepted_move_count: int
    iteration_count: int
    evaluated_neighbor_count: int
    branch_reoptimization_policy: str
    substitution_parameter_policy: str
    substitution_parameter_values: dict[str, float]
    substitution_parameter_warnings: list[str]
    total_branch_optimization_pass_count: int
    total_branch_function_evaluation_count: int
    stopping_reason: str
    equal_best_tree_report: NucleotideLikelihoodEqualBestTreeReport
    trace_rows: list[NucleotideLikelihoodNniTraceRow]
    candidate_rows: list[NucleotideLikelihoodNniCandidateRow]


@dataclass(slots=True)
class NucleotideLikelihoodSprSearchBudget:
    """Normalized rooted likelihood SPR search budget settings."""

    max_candidate_count: int | None
    max_iteration_count: int | None
    max_elapsed_seconds: float | None
    max_accepted_move_count: int | None


@dataclass(slots=True)
class NucleotideLikelihoodSprTraceRow:
    """One deterministic event row in a rooted nucleotide likelihood SPR search trace."""

    event_index: int
    event_kind: str
    iteration: int
    move_type: str
    candidate_topology_fingerprint: str
    log_likelihood_before: float | None
    log_likelihood_after: float
    log_likelihood_delta: float | None
    accepted_move: bool
    trace_reason: str
    tree_before_newick: str | None
    tree_after_newick: str
    pruned_clade_id: str | None
    regraft_target_branch_id: str | None
    branch_reoptimization_policy: str
    branch_reoptimization_scope: str
    affected_branch_clade_ids: list[str]
    optimized_branch_count: int
    optimized_branch_clade_ids: list[str]
    branch_reoptimization_converged: bool | None
    branch_optimization_pass_count: int
    branch_function_evaluation_count: int
    boundary_warning_messages: list[str]
    stopping_reason: str | None
    unsearched_candidate_count: int | None


@dataclass(slots=True)
class NucleotideLikelihoodSprSearchReport:
    """Complete rooted nucleotide likelihood SPR hill-climb report."""

    algorithm: str
    model_name: str
    tree_path: str | None
    alignment_path: str | None
    taxon_count: int
    site_count: int
    pattern_count: int
    input_tree_newick: str
    start_tree_newick: str
    start_log_likelihood: float
    final_tree_newick: str
    final_log_likelihood: float
    accepted_move_count: int
    iteration_count: int
    evaluated_neighbor_count: int
    evaluation_budget: int | None
    search_budget: NucleotideLikelihoodSprSearchBudget
    unsearched_candidate_count: int
    branch_reoptimization_policy: str
    substitution_parameter_policy: str
    substitution_parameter_values: dict[str, float]
    substitution_parameter_warnings: list[str]
    total_branch_optimization_pass_count: int
    total_branch_function_evaluation_count: int
    stopping_reason: str
    equal_best_tree_report: NucleotideLikelihoodEqualBestTreeReport
    trace_rows: list[NucleotideLikelihoodSprTraceRow]


@dataclass(slots=True)
class NucleotideLikelihoodTbrTraceRow:
    """One deterministic event row in a rooted nucleotide likelihood TBR search trace."""

    event_index: int
    event_kind: str
    iteration: int
    move_type: str
    candidate_topology_fingerprint: str
    log_likelihood_before: float | None
    log_likelihood_after: float
    log_likelihood_delta: float | None
    accepted_move: bool
    trace_reason: str
    tree_before_newick: str | None
    tree_after_newick: str
    cut_edge_id: str | None
    left_attachment_branch_id: str | None
    right_attachment_branch_id: str | None
    branch_reoptimization_policy: str
    branch_reoptimization_scope: str
    optimized_branch_count: int
    optimized_branch_clade_ids: list[str]
    branch_reoptimization_converged: bool | None
    branch_optimization_pass_count: int
    branch_function_evaluation_count: int
    boundary_warning_messages: list[str]
    stopping_reason: str | None


@dataclass(slots=True)
class NucleotideLikelihoodTbrSearchReport:
    """Complete rooted nucleotide likelihood TBR hill-climb report."""

    algorithm: str
    model_name: str
    tree_path: str | None
    alignment_path: str | None
    taxon_count: int
    site_count: int
    pattern_count: int
    input_tree_newick: str
    start_tree_newick: str
    start_log_likelihood: float
    final_tree_newick: str
    final_log_likelihood: float
    accepted_move_count: int
    evaluated_neighbor_count: int
    branch_reoptimization_policy: str
    substitution_parameter_policy: str
    substitution_parameter_values: dict[str, float]
    substitution_parameter_warnings: list[str]
    total_branch_optimization_pass_count: int
    total_branch_function_evaluation_count: int
    stopping_reason: str
    equal_best_tree_report: NucleotideLikelihoodEqualBestTreeReport
    trace_rows: list[NucleotideLikelihoodTbrTraceRow]


@dataclass(slots=True)
class NucleotideLikelihoodTopologyPerturbationStep:
    """One random legal topology perturbation applied before native local search."""

    step_index: int
    move_family: str
    tree_before_newick: str
    tree_after_newick: str
    topology_fingerprint_before: str
    topology_fingerprint_after: str
    pivot_branch_id: str | None
    sibling_clade_id: str | None
    exchanged_clade_id: str | None
    pruned_clade_id: str | None
    regraft_target_branch_id: str | None


@dataclass(slots=True)
class NucleotideLikelihoodStochasticTopologyPerturbationSearchReport:
    """Complete stochastic topology-perturbation workflow over one local search run."""

    algorithm: str
    model_name: str
    perturbation_move_family: str
    local_search_method: str
    branch_reoptimization_policy: str
    tree_path: str | None
    alignment_path: str | None
    taxon_count: int
    site_count: int
    pattern_count: int
    perturbation_seed: int
    perturbation_move_count_requested: int
    perturbation_move_count_applied: int
    input_tree_newick: str
    perturbed_tree_newick: str
    perturbed_topology_fingerprint: str
    local_search_algorithm: str
    local_search_start_tree_newick: str
    local_search_start_log_likelihood: float
    final_tree_newick: str
    final_log_likelihood: float
    final_topology_fingerprint: str
    local_search_accepted_move_count: int
    local_search_evaluated_neighbor_count: int
    local_search_stopping_reason: str
    substitution_parameter_policy: str
    substitution_parameter_values: dict[str, float]
    substitution_parameter_warnings: list[str]
    perturbation_steps: list[NucleotideLikelihoodTopologyPerturbationStep]


@dataclass(slots=True)
class NucleotideLikelihoodRatchetCycle:
    """One temporary-reweighting cycle inside a native likelihood ratchet search."""

    cycle_index: int
    start_log_likelihood: float
    start_tree_newick: str
    reweighted_site_positions: list[int]
    temporary_site_weights: dict[int, int]
    perturbation_factor: int
    perturbed_alignment_length: int
    perturbed_pattern_count: int
    perturbed_search_algorithm: str
    perturbed_score: float
    perturbed_tree_newick: str
    perturbed_accepted_move_count: int
    perturbed_evaluated_neighbor_count: int
    perturbed_stopping_reason: str
    restored_search_algorithm: str
    restored_score: float
    restored_tree_newick: str
    restored_accepted_move_count: int
    restored_evaluated_neighbor_count: int
    restored_stopping_reason: str
    best_score_after_cycle: float
    best_tree_after_cycle: str
    best_tree_improved: bool


@dataclass(slots=True)
class NucleotideLikelihoodRatchetBestTreeHistory:
    """One retained best-tree checkpoint inside a native likelihood ratchet search."""

    history_index: int
    cycle_index: int
    best_log_likelihood: float
    best_tree_newick: str
    best_topology_fingerprint: str


@dataclass(slots=True)
class NucleotideLikelihoodRatchetReport:
    """Complete native likelihood ratchet search summary."""

    algorithm: str
    model_name: str
    local_search_method: str
    tree_path: str | None
    alignment_path: str | None
    taxon_count: int
    site_count: int
    pattern_count: int
    cycle_count: int
    perturbation_seed: int
    perturbed_site_count: int
    perturbation_factor: int
    input_tree_newick: str
    start_tree_newick: str
    start_log_likelihood: float
    final_tree_newick: str
    final_log_likelihood: float
    best_tree_newick: str
    best_log_likelihood: float
    branch_reoptimization_policy: str
    evaluation_budget: int | None
    substitution_parameter_policy: str
    substitution_parameter_values: dict[str, float]
    substitution_parameter_warnings: list[str]
    cycle_rows: list[NucleotideLikelihoodRatchetCycle]
    best_tree_history_rows: list[NucleotideLikelihoodRatchetBestTreeHistory]


@dataclass(slots=True)
class NucleotideLikelihoodSimulatedAnnealingTraceRow:
    """One simulated-annealing proposal row inside a native likelihood topology search."""

    iteration: int
    temperature: float
    move_family: str
    current_log_likelihood_before: float
    proposed_log_likelihood: float
    log_likelihood_delta: float
    acceptance_probability: float
    acceptance_uniform_draw: float
    acceptance_decision: str
    accepted_move: bool
    best_tree_improved: bool
    current_tree_before_newick: str
    proposed_tree_newick: str
    current_tree_after_newick: str
    pivot_branch_id: str | None
    sibling_clade_id: str | None
    exchanged_clade_id: str | None
    pruned_clade_id: str | None
    regraft_target_branch_id: str | None
    branch_reoptimization_policy: str
    branch_reoptimization_scope: str
    optimized_branch_count: int
    optimized_branch_clade_ids: list[str]
    branch_reoptimization_converged: bool
    branch_optimization_pass_count: int
    branch_function_evaluation_count: int
    boundary_warning_messages: list[str]


@dataclass(slots=True)
class NucleotideLikelihoodSimulatedAnnealingSearchReport:
    """Complete native likelihood simulated-annealing topology search summary."""

    algorithm: str
    model_name: str
    proposal_move_family: str
    branch_reoptimization_policy: str
    tree_path: str | None
    alignment_path: str | None
    taxon_count: int
    site_count: int
    pattern_count: int
    annealing_seed: int
    iteration_count_requested: int
    iteration_count_completed: int
    initial_temperature: float
    cooling_rate: float
    input_tree_newick: str
    start_tree_newick: str
    start_log_likelihood: float
    final_tree_newick: str
    final_log_likelihood: float
    best_tree_newick: str
    best_log_likelihood: float
    best_topology_fingerprint: str
    accepted_move_count: int
    rejected_move_count: int
    accepted_worse_move_count: int
    substitution_parameter_policy: str
    substitution_parameter_values: dict[str, float]
    substitution_parameter_warnings: list[str]
    total_branch_optimization_pass_count: int
    total_branch_function_evaluation_count: int
    stopping_reason: str
    trace_rows: list[NucleotideLikelihoodSimulatedAnnealingTraceRow]


@dataclass(slots=True)
class NucleotideLikelihoodMultiStartRunSummary:
    """One independently searched start tree inside a multi-start likelihood workflow."""

    search_algorithm: str
    start_tree_source_kind: str
    start_tree_source_label: str
    start_tree_generation_seed: int | None
    start_tree_newick: str
    start_log_likelihood: float
    final_tree_newick: str
    final_log_likelihood: float
    final_topology_fingerprint: str
    search_iteration_count: int
    accepted_move_count: int
    evaluated_neighbor_count: int
    final_likelihood_rank: int
    branch_reoptimization_policy: str
    substitution_parameter_policy: str
    substitution_parameter_values: dict[str, float]
    substitution_parameter_warnings: list[str]
    total_branch_optimization_pass_count: int
    total_branch_function_evaluation_count: int
    stopping_reason: str
    best_run: bool


@dataclass(slots=True)
class NucleotideLikelihoodMultiStartSearchReport:
    """Complete rooted nucleotide likelihood multi-start search summary."""

    algorithm: str
    model_name: str
    local_search_method: str
    tree_path: str | None
    alignment_path: str | None
    taxon_count: int
    site_count: int
    pattern_count: int
    input_tree_newick: str
    start_tree_source_policy: str
    starting_tree_selection_policy: str | None
    input_tree_included: bool
    available_start_tree_count: int
    generated_start_tree_count: int
    start_tree_count: int
    start_tree_seed: int
    evaluation_budget: int | None
    branch_reoptimization_policy: str
    best_run_source_label: str
    best_final_tree_newick: str
    best_final_log_likelihood: float
    best_final_topology_fingerprint: str
    run_summaries: list[NucleotideLikelihoodMultiStartRunSummary]


@dataclass(slots=True)
class NucleotideLikelihoodTreeInferenceReport:
    """Complete native nucleotide maximum-likelihood tree inference summary."""

    algorithm: str
    alignment_path: str | None
    taxon_count: int
    site_count: int
    pattern_count: int
    stepwise_addition_model_name: str
    stepwise_addition_tree_newick: str
    stepwise_addition_final_score: float
    start_tree_source_policy: str
    random_start_tree_count: int
    start_tree_seed: int
    model_selection_strategy: str
    model_selection_criterion: str | None
    model_selection_tree_newick: str
    selected_model_name: str
    search_method: str
    branch_reoptimization_policy: str
    run_summaries: list[NucleotideLikelihoodMultiStartRunSummary]
    best_run_source_label: str
    best_final_tree_newick: str
    best_final_log_likelihood: float
    best_final_topology_fingerprint: str
    model_selection_report: SubstitutionModelSelectionReport
    best_search_report: (
        NucleotideLikelihoodNniSearchReport
        | NucleotideLikelihoodSprSearchReport
        | NucleotideLikelihoodTbrSearchReport
    )


@dataclass(slots=True)
class LikelihoodWrapperCorrespondenceObservation:
    """One governed native-versus-wrapper ML correspondence case."""

    case_id: str
    wrapper_engine: str
    native_surface: str
    wrapper_surface: str
    comparison_policy: str
    status: str
    supported: bool
    blocking: bool
    tolerance: float | None
    rationale: str
    input_fixtures: list[str]
    expected_output: dict[str, object]
    observed_output: dict[str, object]


@dataclass(slots=True)
class LikelihoodWrapperCorrespondenceSummaryRow:
    """One status-level summary row across governed ML correspondence cases."""

    status: str
    case_count: int
    blocking_case_count: int
    case_ids: list[str]


@dataclass(slots=True)
class LikelihoodWrapperCorrespondenceReport:
    """One governed correspondence report across cached ML wrapper references."""

    observations: list[LikelihoodWrapperCorrespondenceObservation]
    summary_rows: list[LikelihoodWrapperCorrespondenceSummaryRow]
    case_count: int
    supported_case_count: int
    exact_match_case_count: int
    tolerance_match_case_count: int
    expected_model_assumption_difference_case_count: int
    unsupported_case_count: int
    native_bug_case_count: int
    blocking_case_count: int
    all_supported_cases_clear: bool


@dataclass(slots=True)
class NucleotideLikelihoodBootstrapReplicateRow:
    """One site-resampled native ML tree-inference replicate."""

    replicate_index: int
    sampled_site_indices: list[int]
    replicate_start_tree_seed: int
    selected_model_name: str
    best_run_source_label: str
    final_tree_newick: str
    final_log_likelihood: float
    final_topology_fingerprint: str
    accepted_move_count: int
    search_iteration_count: int


@dataclass(slots=True)
class NucleotideLikelihoodBootstrapCladeSupportRow:
    """One reference-tree bootstrap support row mapped by descendant-tip set."""

    branch_id: str
    node_label: str | None
    descendant_taxa: list[str]
    supporting_tree_count: int
    clade_frequency: float
    support_percent: float


@dataclass(slots=True)
class NucleotideLikelihoodBootstrapTreeInferenceReport:
    """Native nucleotide ML bootstrap inference summary over one alignment."""

    algorithm: str
    alignment_path: str | None
    requested_model_name: str
    model_selection_strategy: str
    model_selection_criterion: str | None
    selected_reference_model_name: str
    search_method: str
    branch_reoptimization_policy: str
    taxon_count: int
    site_count: int
    pattern_count: int
    start_tree_count: int
    start_tree_seed: int
    replicate_count: int
    bootstrap_seed: int
    reference_tree_newick: str
    reference_log_likelihood: float
    reference_topology_fingerprint: str
    reference_best_run_source_label: str
    replicate_rows: list[NucleotideLikelihoodBootstrapReplicateRow]
    clade_support_rows: list[NucleotideLikelihoodBootstrapCladeSupportRow]


@dataclass(slots=True)
class NucleotideShLikeBranchSupportRow:
    """One native SH-like local support summary for a rooted internal branch."""

    branch_id: str
    node_label: str | None
    descendant_taxa: list[str]
    alternative_arrangement_count: int
    reference_log_likelihood: float
    best_alternative_tree_id: str
    best_alternative_tree_label: str
    best_alternative_topology_fingerprint: str
    best_alternative_log_likelihood: float
    observed_delta_log_likelihood: float
    reference_is_observed_best: bool
    support_replicate_count: int
    support_fraction: float
    support_percent: float
    caution_label: str


@dataclass(slots=True)
class NucleotideShLikeBranchSupportLocalTopologyRow:
    """One local rooted-NNI arrangement scored for SH-like branch support."""

    branch_id: str
    node_label: str | None
    descendant_taxa: list[str]
    candidate_tree_id: str
    candidate_tree_label: str
    local_arrangement_kind: str
    tree_newick: str
    topology_fingerprint: str
    observed_log_likelihood: float
    observed_delta_log_likelihood: float
    observed_best_local_arrangement: bool


@dataclass(slots=True)
class NucleotideShLikeBranchSupportResamplingRow:
    """One site-resampled local branch-support comparison against NNI alternatives."""

    branch_id: str
    descendant_taxa: list[str]
    replicate_index: int
    reference_resampled_log_likelihood: float
    best_local_tree_id: str
    best_local_tree_label: str
    best_local_resampled_log_likelihood: float
    best_alternative_tree_id: str
    best_alternative_tree_label: str
    best_alternative_resampled_log_likelihood: float
    reference_delta_log_likelihood: float
    reference_matches_or_beats_alternatives: bool


@dataclass(slots=True)
class NucleotideShLikeBranchSupportReport:
    """Native SH-like approximate local branch support over one rooted tree."""

    algorithm: str
    tree_path: str | None
    alignment_path: str | None
    model_name: str
    taxa: list[str]
    branch_count: int
    site_count: int
    pattern_count: int
    compression_used: bool
    parameter_values: dict[str, float]
    reference_tree_newick: str
    resampling_method: str
    resampling_replicate_count: int
    resampling_seed: int
    caution_label: str
    branch_support_rows: list[NucleotideShLikeBranchSupportRow]
    local_topology_rows: list[NucleotideShLikeBranchSupportLocalTopologyRow]
    resampling_rows: list[NucleotideShLikeBranchSupportResamplingRow]


@dataclass(slots=True)
class NucleotideLikelihoodSearchTraceReplayStep:
    """One accepted search-trace step replayed against the evolving tree."""

    trace_event_index: int
    iteration: int
    move_type: str
    trace_reason: str
    candidate_topology_fingerprint: str
    replayed_topology_fingerprint: str | None
    matched_candidate_count: int
    step_replayed: bool


@dataclass(slots=True)
class NucleotideLikelihoodSearchTraceReplayReport:
    """Replay summary for one native likelihood local-search trace."""

    algorithm: str
    source_search_algorithm: str
    accepted_trace_event_count: int
    replayed_step_count: int
    start_tree_newick: str
    stored_final_tree_newick: str
    stored_final_topology_fingerprint: str
    replayed_final_tree_newick: str | None
    replayed_final_topology_fingerprint: str | None
    final_topology_matches: bool
    replay_failed: bool
    failure_reason: str | None
    step_rows: list[NucleotideLikelihoodSearchTraceReplayStep]


@dataclass(slots=True)
class NucleotideLikelihoodStartingTreeSummary:
    """One scored starting tree inside a native likelihood start-tree pool."""

    tree_id: str
    source_strategy: str
    generation_seed: int | None
    topology_hash: str
    starting_log_likelihood: float
    substitution_parameter_policy: str
    substitution_parameter_values: dict[str, float]
    substitution_parameter_warnings: list[str]
    tree_newick: str


@dataclass(slots=True)
class NucleotideLikelihoodStartingTreePoolReport:
    """Scored pool of distinct likelihood starting trees for one alignment."""

    algorithm: str
    model_name: str
    tree_path: str | None
    alignment_path: str | None
    taxon_count: int
    site_count: int
    pattern_count: int
    random_start_tree_count: int
    random_start_tree_seed: int
    starting_tree_summaries: list[NucleotideLikelihoodStartingTreeSummary]


@dataclass(slots=True)
class CandidateTreeSiteLikelihoodSummary:
    """One candidate-tree summary inside a shared site-likelihood matrix workflow."""

    candidate_tree_id: str
    candidate_tree_label: str
    topology_fingerprint: str
    tree_newick: str
    log_likelihood: float
    observed_delta_log_likelihood: float


@dataclass(slots=True)
class CandidateTreeSiteLikelihoodRow:
    """One candidate-tree by site log-likelihood matrix row."""

    candidate_tree_id: str
    candidate_tree_label: str
    tree_newick: str
    pattern_id: str
    pattern_weight: int
    site_position: int
    site_states: tuple[str, ...]
    log_likelihood: float


@dataclass(slots=True)
class CandidateTreeSiteLikelihoodMatrixReport:
    """Expanded candidate-tree by site likelihood matrix for one shared alignment/model."""

    model_name: str
    tree_set_path: str | None
    alignment_path: str | None
    taxa: list[str]
    tree_count: int
    site_count: int
    pattern_count: int
    compression_used: bool
    expansion_policy: str
    comparison_caution_label: str
    parameter_values: dict[str, float]
    candidate_trees: list[CandidateTreeSiteLikelihoodSummary]
    matrix_rows: list[CandidateTreeSiteLikelihoodRow]


@dataclass(slots=True)
class ApproximateTopologyTestSummaryRow:
    """One candidate-tree summary row inside a site-resampled topology comparison."""

    candidate_tree_id: str
    candidate_tree_label: str
    tree_newick: str
    observed_log_likelihood: float
    observed_delta_log_likelihood: float
    observed_best_tree: bool
    resampling_win_count: int
    resampling_frequency: float
    p_like_statistic: float
    resampling_mean_delta_log_likelihood: float
    resampling_min_delta_log_likelihood: float
    resampling_max_delta_log_likelihood: float
    caution_label: str


@dataclass(slots=True)
class ApproximateTopologyTestResamplingRow:
    """One resampled candidate-tree delta row in an approximate topology comparison."""

    replicate_index: int
    candidate_tree_id: str
    candidate_tree_label: str
    resampled_log_likelihood: float
    observed_best_tree_id: str
    observed_best_tree_label: str
    observed_best_resampled_log_likelihood: float
    resampled_delta_log_likelihood: float
    candidate_matches_or_beats_observed_best: bool


@dataclass(slots=True)
class ApproximateTopologyTestReport:
    """Site-resampled approximate topology comparison across candidate trees."""

    algorithm: str
    model_name: str
    tree_set_path: str | None
    alignment_path: str | None
    taxa: list[str]
    tree_count: int
    site_count: int
    pattern_count: int
    compression_used: bool
    expansion_policy: str
    parameter_values: dict[str, float]
    resampling_method: str
    resampling_replicate_count: int
    resampling_seed: int
    observed_best_tree_id: str
    observed_best_tree_label: str
    caution_label: str
    summary_rows: list[ApproximateTopologyTestSummaryRow]
    resampling_rows: list[ApproximateTopologyTestResamplingRow]


@dataclass(slots=True)
class SiteLogLikelihoodRow:
    """One expanded site log-likelihood row from one compressed-pattern evaluation."""

    pattern_id: str
    pattern_weight: int
    site_position: int
    site_states: tuple[str, ...]
    log_likelihood: float


@dataclass(slots=True)
class FixedTopologySiteLogLikelihoodReport:
    """Per-site fixed-topology likelihood report with explicit expansion policy."""

    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    expansion_policy: str
    tree_newick: str
    state_count: int
    observation_policy: str
    parameter_values: dict[str, float]
    log_likelihood: float
    site_log_likelihoods: list[SiteLogLikelihoodRow]


@dataclass(slots=True)
class BranchLikelihoodDiagnosticRow:
    """One branchwise likelihood diagnostic derived from branch-collapse replay."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    branch_length: float
    collapsed_branch_log_likelihood: float
    contribution_proxy: float
    warning_flags: list[str]


@dataclass(slots=True)
class FixedTreeBranchLikelihoodDiagnosticsReport:
    """Fixed-tree branchwise likelihood diagnostics under one sequence model."""

    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    branch_count: int
    tree_newick: str
    baseline_log_likelihood: float
    branch_diagnostics: list[BranchLikelihoodDiagnosticRow]


@dataclass(slots=True)
class MarginalAncestralStateProbabilityRow:
    """One internal-node posterior probability for one site and one state."""

    node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    pattern_id: str
    site_position: int
    state: str
    posterior_probability: float


@dataclass(slots=True)
class MarginalAncestralSequenceProbabilityReport:
    """Expanded internal-node marginal sequence posterior report."""

    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    internal_node_count: int
    compression_used: bool
    expansion_policy: str
    tree_newick: str
    parameter_values: dict[str, float]
    posterior_rows: list[MarginalAncestralStateProbabilityRow]


@dataclass(slots=True)
class MarginalAncestralSiteSummaryRow:
    """One internal-node posterior summary for one site across all DNA states."""

    node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    pattern_id: str
    site_position: int
    most_likely_state: str
    max_posterior_probability: float
    posterior_probability_a: float
    posterior_probability_c: float
    posterior_probability_g: float
    posterior_probability_t: float


@dataclass(slots=True)
class MarginalAncestralSequenceExportRecord:
    """One FASTA-ready ancestral sequence for one internal node."""

    node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    sequence: str


@dataclass(slots=True)
class MarginalAncestralSequenceUncertaintyRow:
    """One node-site export summary under one explicit low-confidence policy."""

    node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    pattern_id: str
    site_position: int
    exported_state: str
    most_likely_state: str
    max_posterior_probability: float
    low_confidence: bool
    posterior_probability_a: float
    posterior_probability_c: float
    posterior_probability_g: float
    posterior_probability_t: float


@dataclass(slots=True)
class MarginalAncestralSequenceFastaExportReport:
    """FASTA-oriented marginal ancestral sequence export with uncertainty rows."""

    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    internal_node_count: int
    compression_used: bool
    expansion_policy: str
    tree_newick: str
    parameter_values: dict[str, float]
    posterior_probability_threshold: float
    low_confidence_state_symbol: str
    sequence_records: list[MarginalAncestralSequenceExportRecord]
    uncertainty_rows: list[MarginalAncestralSequenceUncertaintyRow]


@dataclass(slots=True)
class JointAncestralStateAssignmentRow:
    """One joint ancestral state assignment for one internal node and one site."""

    node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    pattern_id: str
    site_position: int
    state: str


@dataclass(slots=True)
class JointAncestralSequenceRecord:
    """One fixed-topology joint ancestral sequence for one internal node."""

    node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    sequence: str


@dataclass(slots=True)
class JointAncestralSequenceReport:
    """Expanded fixed-topology joint ancestral sequence reconstruction report."""

    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    internal_node_count: int
    compression_used: bool
    expansion_policy: str
    tree_newick: str
    parameter_values: dict[str, float]
    sequence_records: list[JointAncestralSequenceRecord]
    assignment_rows: list[JointAncestralStateAssignmentRow]


@dataclass(slots=True)
class ProteinPoissonTreeLikelihoodReport:
    """Native 20-state protein Poisson likelihood report for one fixed topology."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    state_count: int
    gap_policy: str
    missing_policy: str
    ambiguity_policy: str
    log_likelihood: float
    site_log_likelihoods: list[SiteLogLikelihoodRow]
    branch_likelihood_diagnostics: list[BranchLikelihoodDiagnosticRow] | None = None


@dataclass(slots=True)
class CodonCtmcTreeLikelihoodReport:
    """Native sense-codon CTMC likelihood report for one fixed topology."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    state_count: int
    genetic_code_id: int
    genetic_code_name: str
    codon_frequency_source: str
    log_likelihood: float


@dataclass(slots=True)
class ProteinEmpiricalMatrixTreeLikelihoodReport:
    """Native empirical 20-state protein likelihood report for one fixed topology."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    state_count: int
    matrix_label: str
    root_prior_source: str
    gap_policy: str
    missing_policy: str
    log_likelihood: float
    site_log_likelihoods: list[SiteLogLikelihoodRow]
    branch_likelihood_diagnostics: list[BranchLikelihoodDiagnosticRow] | None = None


@dataclass(slots=True)
class BranchLengthOptimizationRow:
    """One branch-length change recorded from one fixed-topology optimization run."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    initial_branch_length: float
    optimized_branch_length: float


@dataclass(slots=True)
class FixedTopologyNucleotideBranchLengthOptimizationReport:
    """Fixed-topology branch-length optimization summary for one nucleotide model."""

    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    branch_count: int
    initial_tree_newick: str
    optimized_tree_newick: str
    state_count: int
    observation_policy: str
    root_prior_source: str
    root_prior_values: list[float]
    parameter_count: int
    fixed_parameter_values: dict[str, float]
    initial_log_likelihood: float
    optimized_log_likelihood: float
    optimization_pass_count: int
    function_evaluation_count: int
    converged: bool
    lower_branch_length_bound: float
    upper_branch_length_bound: float
    branches: list[BranchLengthOptimizationRow]
    boundary_warnings: list[LikelihoodOptimizationBoundaryWarning]


@dataclass(slots=True)
class FixedTopologyNucleotideSingleBranchOptimizationReport:
    """Fixed-topology optimization summary for one selected nucleotide branch."""

    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    branch_count: int
    initial_tree_newick: str
    optimized_tree_newick: str
    state_count: int
    observation_policy: str
    root_prior_source: str
    root_prior_values: list[float]
    parameter_count: int
    fixed_parameter_values: dict[str, float]
    selected_branch: BranchLengthOptimizationRow
    unchanged_branch_count: int
    unchanged_branch_ids: list[str]
    initial_log_likelihood: float
    optimized_log_likelihood: float
    function_evaluation_count: int
    converged: bool
    lower_branch_length_bound: float
    upper_branch_length_bound: float


@dataclass(slots=True)
class JointNucleotideOptimizationUpdateRow:
    """One branch or substitution update inside a joint fixed-topology optimization run."""

    joint_pass_index: int
    update_kind: str
    starting_log_likelihood: float
    optimized_log_likelihood: float
    log_likelihood_delta: float
    function_evaluation_count: int
    optimization_pass_count: int
    converged: bool
    optimized_branch_count: int
    optimized_branch_ids: list[str]
    updated_parameter_names: list[str]


@dataclass(slots=True)
class FixedTopologyNucleotideJointOptimizationReport:
    """Joint fixed-topology branch-and-model optimization summary for one nucleotide model."""

    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    branch_count: int
    initial_tree_newick: str
    optimized_tree_newick: str
    state_count: int
    observation_policy: str
    root_prior_source: str
    root_prior_values: list[float]
    parameter_count: int
    base_frequency_source: str | None
    base_frequency_a: float | None
    base_frequency_c: float | None
    base_frequency_g: float | None
    base_frequency_t: float | None
    fixed_parameter_values: dict[str, float]
    parameter_rows: list[SubstitutionParameterOptimizationRow]
    branch_rows: list[BranchLengthOptimizationRow]
    initial_log_likelihood: float
    optimized_log_likelihood: float
    function_evaluation_count: int
    joint_optimization_pass_count: int
    converged: bool
    convergence_reason: str
    lower_branch_length_bound: float
    upper_branch_length_bound: float
    update_rows: list[JointNucleotideOptimizationUpdateRow]
    boundary_warnings: list[LikelihoodOptimizationBoundaryWarning]
    warnings: list[str]


@dataclass(slots=True)
class JointNucleotideOptimizationRestartAttemptRow:
    """One recorded restart attempt around one joint branch-and-model optimization run."""

    attempt_index: int
    trigger_reason: str
    initial_kappa: float | None
    initial_exchangeability_profile_name: str | None
    optimized_log_likelihood: float
    converged: bool
    convergence_reason: str
    boundary_warning_count: int
    branch_boundary_count: int
    warning_count: int
    selected_best: bool


@dataclass(slots=True)
class FixedTopologyNucleotideJointOptimizationRestartReport:
    """Restart-policy summary around one joint fixed-topology nucleotide optimization workflow."""

    model_name: str
    restart_policy: str
    attempt_count: int
    selected_attempt_index: int
    selected_solution_reason: str
    selected_report: FixedTopologyNucleotideJointOptimizationReport
    attempt_rows: list[JointNucleotideOptimizationRestartAttemptRow]


@dataclass(slots=True)
class NucleotideLikelihoodOptimizationEquivalenceReport:
    """Independent rescore comparison for one optimized fixed-topology nucleotide result."""

    optimization_surface: str
    model_name: str
    taxa: list[str]
    site_count: int
    pattern_count: int
    optimized_tree_newick: str
    parameter_values: dict[str, float]
    root_prior_source: str | None
    stored_log_likelihood: float
    independently_rescored_log_likelihood: float
    absolute_difference: float
    relative_difference: float
    absolute_tolerance: float
    relative_tolerance: float
    equivalent: bool


@dataclass(slots=True)
class ProteinEmpiricalBranchLengthOptimizationReport:
    """Fixed-topology branch-length optimization summary for one empirical protein model."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    branch_count: int
    initial_tree_newick: str
    optimized_tree_newick: str
    state_count: int
    matrix_label: str
    root_prior_source: str
    gap_policy: str
    missing_policy: str
    likelihood_model: str
    alpha: float | None
    invariant_proportion: float | None
    initial_log_likelihood: float
    optimized_log_likelihood: float
    optimization_pass_count: int
    function_evaluation_count: int
    converged: bool
    lower_branch_length_bound: float
    upper_branch_length_bound: float
    branches: list[BranchLengthOptimizationRow]


@dataclass(slots=True)
class DiscreteGammaRateCategory:
    """One discrete-gamma rate category with one stable weight."""

    category_index: int
    rate: float
    weight: float


@dataclass(slots=True)
class DiscreteGammaSiteLikelihood:
    """One emitted site likelihood row under one gamma-mixture evaluation."""

    pattern_id: str
    pattern_weight: int
    site_position: int
    site_states: tuple[str, ...]
    category_likelihoods: list[float]
    mixture_likelihood: float
    log_likelihood: float


@dataclass(slots=True)
class ProteinEmpiricalDiscreteGammaTreeLikelihoodReport:
    """Empirical protein likelihood report with discrete-gamma rate heterogeneity."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    state_count: int
    matrix_label: str
    root_prior_source: str
    gap_policy: str
    missing_policy: str
    alpha: float
    category_count: int
    category_rates: list[DiscreteGammaRateCategory]
    site_likelihoods: list[DiscreteGammaSiteLikelihood]
    log_likelihood: float
    branch_likelihood_diagnostics: list[BranchLikelihoodDiagnosticRow] | None = None


@dataclass(slots=True)
class InvariantMixtureSiteLikelihood:
    """One emitted site likelihood row under one invariant-site mixture evaluation."""

    pattern_id: str
    pattern_weight: int
    site_position: int
    site_states: tuple[str, ...]
    invariant_component_likelihood: float
    variable_component_likelihood: float
    mixture_likelihood: float
    log_likelihood: float


@dataclass(slots=True)
class ProteinEmpiricalInvariantMixtureTreeLikelihoodReport:
    """Empirical protein likelihood report with one fitted invariant-site mixture."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    state_count: int
    matrix_label: str
    root_prior_source: str
    gap_policy: str
    missing_policy: str
    initial_invariant_proportion: float
    invariant_proportion: float
    initial_log_likelihood: float
    log_likelihood: float
    function_evaluation_count: int
    converged: bool
    lower_invariant_proportion_bound: float
    upper_invariant_proportion_bound: float
    site_likelihoods: list[InvariantMixtureSiteLikelihood]
    branch_likelihood_diagnostics: list[BranchLikelihoodDiagnosticRow] | None = None


@dataclass(slots=True)
class DiscreteGammaInvariantMixtureSiteLikelihood:
    """One emitted site row under one combined discrete-gamma invariant mixture."""

    pattern_id: str
    pattern_weight: int
    site_position: int
    site_states: tuple[str, ...]
    category_likelihoods: list[float]
    invariant_component_likelihood: float
    variable_component_likelihood: float
    mixture_likelihood: float
    log_likelihood: float


@dataclass(slots=True)
class ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport:
    """Empirical protein likelihood report with active gamma and invariant mixture."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
    state_count: int
    matrix_label: str
    root_prior_source: str
    gap_policy: str
    missing_policy: str
    alpha: float
    category_count: int
    category_rates: list[DiscreteGammaRateCategory]
    initial_invariant_proportion: float
    invariant_proportion: float
    initial_log_likelihood: float
    log_likelihood: float
    function_evaluation_count: int
    converged: bool
    lower_invariant_proportion_bound: float
    upper_invariant_proportion_bound: float
    hit_lower_invariant_proportion_boundary: bool
    hit_upper_invariant_proportion_boundary: bool
    boundary_warnings: list[str]
    site_likelihoods: list[DiscreteGammaInvariantMixtureSiteLikelihood]
    branch_likelihood_diagnostics: list[BranchLikelihoodDiagnosticRow] | None = None


@dataclass(frozen=True, slots=True)
class LikelihoodPlacementAlternativeRow:
    """One optimized query placement on one candidate reference-tree edge."""

    query_id: str
    placement_rank: int
    edge_id: str
    child_name: str | None
    descendant_taxa: list[str]
    original_branch_length: float
    optimized_proximal_length: float
    optimized_distal_length: float
    optimized_pendant_length: float
    log_likelihood: float
    likelihood_weight_ratio: float
    function_evaluation_count: int
    optimization_pass_count: int
    converged: bool
    placed_tree_newick: str


@dataclass(frozen=True, slots=True)
class LikelihoodPlacementQuerySummary:
    """Best-placement summary for one query sequence."""

    query_id: str
    site_count: int
    pattern_count: int
    best_edge_id: str
    best_child_name: str | None
    best_descendant_taxa: list[str]
    best_original_branch_length: float
    best_proximal_length: float
    best_distal_length: float
    best_pendant_length: float
    best_log_likelihood: float
    best_likelihood_weight_ratio: float
    candidate_placement_count: int
    equally_best_placement_count: int
    best_tree_newick: str


@dataclass(slots=True)
class LikelihoodPlacementReport:
    """Likelihood placement summary for multiple query sequences on one reference tree."""

    model_name: str
    reference_tree_newick: str
    reference_taxa: list[str]
    edge_count: int
    query_count: int
    site_count: int
    lower_pendant_length_bound: float
    upper_pendant_length_bound: float
    max_coordinate_passes: int
    total_function_evaluation_count: int
    query_summaries: list[LikelihoodPlacementQuerySummary]
    alternative_placements: list[LikelihoodPlacementAlternativeRow]
