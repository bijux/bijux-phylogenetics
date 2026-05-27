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


@dataclass(slots=True)
class K80TreeLikelihoodReport:
    """Native K80 likelihood report for one fixed topology and alignment."""

    taxa: list[str]
    site_count: int
    pattern_count: int
    compression_used: bool
    tree_newick: str
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
    warnings: list[str]


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
    parameter_values: dict[str, float]
    log_likelihood: float
    site_log_likelihoods: list[SiteLogLikelihoodRow]


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


@dataclass(slots=True)
class BranchLengthOptimizationRow:
    """One branch-length change recorded from one fixed-topology optimization run."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    initial_branch_length: float
    optimized_branch_length: float


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
    site_position: int
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


@dataclass(slots=True)
class InvariantMixtureSiteLikelihood:
    """One emitted site likelihood row under one invariant-site mixture evaluation."""

    pattern_id: str
    site_position: int
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


@dataclass(slots=True)
class DiscreteGammaInvariantMixtureSiteLikelihood:
    """One emitted site row under one combined discrete-gamma invariant mixture."""

    pattern_id: str
    site_position: int
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
