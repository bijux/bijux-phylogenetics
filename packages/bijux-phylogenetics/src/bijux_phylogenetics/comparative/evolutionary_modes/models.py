from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.continuous.model_fitting import (
    ComparativeResidualSummary,
)
from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonRow,
)

ALLOWED_EVOLUTIONARY_MODES = {
    "brownian",
    "white-noise",
    "pagel-lambda",
    "pagel-kappa",
    "pagel-delta",
    "ornstein-uhlenbeck",
    "early-burst",
}

EXCLUDED_GEIGER_TREND_MODE_ALIASES = {
    "trend",
    "mean-trend",
    "mean_trend",
    "rate-trend",
    "rate_trend",
}

EXCLUDED_GEIGER_STANDARD_ERROR_POLICY = (
    "fitcontinuous-standard-error-explicitly-excluded-this-round"
)

FITCONTINUOUS_MODEL_COMPARISON_ORDER = (
    "brownian",
    "white-noise",
    "pagel-lambda",
    "pagel-kappa",
    "pagel-delta",
    "ornstein-uhlenbeck",
    "early-burst",
)

CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY = (
    "full-gaussian-loglikelihood-includes-normalizing-constant"
)
CONTINUOUS_GAUSSIAN_LIKELIHOOD_COMPARISON_POLICY = "raw-loglikelihood-and-derived-aic-are-directly-comparable-when-the-shared-gaussian-constant-policy-matches"
FITCONTINUOUS_MODEL_RANKING_POLICY = "relative-aic-and-aicc-ranking-is-permitted-only-when-all-candidate-modes-share-one-gaussian-likelihood-constant-policy"
FITCONTINUOUS_MODEL_CONFIDENCE_WEIGHT_BASIS = "AICc"
FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD = 2.0


@dataclass(slots=True)
class EvolutionaryModeBranchLengthRow:
    """One deterministic branch-length change under a governed tree rescaling."""

    node: str
    descendant_taxa: list[str]
    original_branch_length: float
    transformed_branch_length: float
    parent_depth: float
    child_depth: float


@dataclass(slots=True)
class ComparativeTreeRescalingReport:
    """Canonical summary of one shared geiger-style tree-rescaling surface."""

    tree_path: Path
    mode: str
    parameter_name: str
    parameter_value: float
    tip_count: int
    original_total_branch_length: float
    transformed_total_branch_length: float
    transformed_tree_newick: str
    branch_rows: list[EvolutionaryModeBranchLengthRow]


@dataclass(slots=True)
class EvolutionaryModeIdentifiabilityWarning:
    """One identifiability or boundary warning on an evolutionary-mode fit."""

    kind: str
    message: str


@dataclass(slots=True)
class ContinuousModeOptimizerDiagnostics:
    """Optimizer diagnostics for one governed evolutionary-mode parameter search."""

    optimizer_name: str
    parameter_search_strategy: str
    lower_bound: float
    upper_bound: float
    starting_parameter_policy: str
    starting_parameter_value: float
    starting_parameter_log_likelihood: float
    coarse_grid_point_count: int
    fine_grid_point_count: int
    refinement_start_count: int
    function_evaluation_count: int
    coarse_best_parameter: float
    coarse_best_log_likelihood: float
    fine_search_start: float
    fine_search_stop: float
    converged: bool
    hit_lower_boundary: bool
    hit_upper_boundary: bool


@dataclass(slots=True)
class ContinuousModeOptimizerProfileRow:
    """One evaluated point on a governed single-parameter optimizer surface."""

    parameter_value: float
    log_likelihood: float


@dataclass(frozen=True, slots=True)
class ContinuousModeSearchControls:
    """User-visible bounded-search controls for parameterized continuous-mode fits."""

    coarse_grid_point_count: int = 81
    fine_grid_point_count: int = 81
    initial_parameter_value: float | None = None
    refinement_start_count: int = 3


@dataclass(slots=True)
class ContinuousModeBoundaryAssessment:
    """Boundary-focused interpretation surface for one parameterized fitContinuous fit."""

    affected_parameter: str
    parameter_value: float
    lower_bound: float
    upper_bound: float
    hit_lower_boundary: bool
    hit_upper_boundary: bool
    near_lower_boundary: bool
    near_upper_boundary: bool
    flat_likelihood_near_boundary: bool
    boundary_warning_kinds: list[str]
    boundary_dominates_interpretation: bool
    stable_conclusion_supported: bool


@dataclass(slots=True)
class ContinuousEvolutionaryModeFitReport:
    """Intercept-only continuous-trait fit under one governed evolutionary mode."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    taxon_count: int
    taxa: list[str]
    mode: str
    parameter_name: str | None
    parameter_value: float | None
    root_state: float
    rate: float
    log_likelihood: float
    aic: float
    aicc: float
    likelihood_constant_policy: str
    likelihood_comparison_policy: str
    fitted_values: list[float]
    residuals: list[float]
    transformed_tree_newick: str
    confidence_intervals: list[object]
    residual_diagnostics: ComparativeResidualSummary
    optimizer_diagnostics: ContinuousModeOptimizerDiagnostics | None
    optimizer_profile_rows: list[ContinuousModeOptimizerProfileRow] | None
    identifiability_warnings: list[EvolutionaryModeIdentifiabilityWarning]
    assumptions: list[str]
    boundary_assessment: ContinuousModeBoundaryAssessment | None = None


@dataclass(slots=True)
class LikelihoodRatioTestResult:
    """Likelihood-ratio comparison between two governed evolutionary-mode fits."""

    comparison_id: str
    left_mode: str
    right_mode: str
    statistic: float
    degrees_of_freedom: int
    p_value: float


@dataclass(slots=True)
class ContinuousEvolutionaryModeComparisonReport:
    """Model-comparison summary over governed `fitContinuous`-style mode fits."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    rows: list[ComparativeModelComparisonRow]
    better_model: str
    likelihood_constant_policy: str | None
    likelihood_comparison_policy: str
    noncomparable_likelihood_models: list[str]
    likelihood_ratio_tests: list[LikelihoodRatioTestResult]
    model_confidence_weight_basis: str
    model_confidence_delta_threshold: float
    selected_model_akaike_weight: float | None
    models_within_delta_aic_threshold: list[str]
    models_within_delta_aicc_threshold: list[str]
    uncertainty_language: str
    warnings: list[str]
    selected_model_boundary_assessment: ContinuousModeBoundaryAssessment | None = None
    stable_conclusion_supported: bool = True
