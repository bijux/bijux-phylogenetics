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
