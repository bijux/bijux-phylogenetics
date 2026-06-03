from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import (
    DiscreteModelBaselineComparison,
    DiscreteOptimizerDiagnostics,
    DiscreteTransitionRateRow,
)
from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonRow,
)

DISCRETE_MK_MODEL_COMPARISON_ORDER = (
    "equal-rates",
    "symmetric",
    "all-rates-different",
)
DISCRETE_MK_ASCERTAINMENT_POLICIES = (
    "none",
    "lewis-variable-only",
)
DISCRETE_MK_LIKELIHOOD_CONSTANT_POLICY = (
    "continuous-time-markov-pruning-loglikelihood-has-no-extra-normalizing-constant"
)
DISCRETE_MK_LIKELIHOOD_VARIABLE_ONLY_CONSTANT_POLICY = "continuous-time-markov-pruning-loglikelihood-conditions-on-variable-site-observation-under-the-declared-lewis-mk-ascertainment-policy"
DISCRETE_MK_LIKELIHOOD_COMPARISON_POLICY = "raw-loglikelihood-and-derived-aic-are-directly-comparable-when-all-candidate-mk-models-share-the-owned-pruning-likelihood-policy"
DISCRETE_MK_MODEL_RANKING_POLICY = "relative-aic-and-aicc-ranking-is-permitted-only-when-all-candidate-discrete-mk-models-share-one-pruning-likelihood-policy"
DISCRETE_MK_MODEL_CONFIDENCE_WEIGHT_BASIS = "AICc"
DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD = 2.0


def validate_discrete_mk_ascertainment_policy(policy: str) -> str:
    """Normalize and validate one owned discrete Mk ascertainment policy."""
    normalized_policy = policy.strip().lower()
    if normalized_policy not in DISCRETE_MK_ASCERTAINMENT_POLICIES:
        expected = ", ".join(DISCRETE_MK_ASCERTAINMENT_POLICIES)
        raise ValueError(
            "unsupported discrete Mk ascertainment policy: "
            f"{policy}; expected one of {expected}"
        )
    return normalized_policy


def resolve_discrete_mk_likelihood_constant_policy(
    ascertainment_policy: str,
) -> str:
    """Return the governed likelihood-constant policy for one Mk fit surface."""
    validated_policy = validate_discrete_mk_ascertainment_policy(ascertainment_policy)
    if validated_policy == "lewis-variable-only":
        return DISCRETE_MK_LIKELIHOOD_VARIABLE_ONLY_CONSTANT_POLICY
    return DISCRETE_MK_LIKELIHOOD_CONSTANT_POLICY


@dataclass(slots=True)
class DiscreteMkInputAudit:
    """Owned input-policy audit for one discrete Mk model fit."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    taxa: list[str]
    observed_states: list[str]
    state_counts: dict[str, int]
    sparse_states: list[str]
    tree_is_ultrametric: bool
    minimum_root_to_tip_depth: float
    maximum_root_to_tip_depth: float
    ultrametric_policy: str
    missing_value_policy: str
    missing_from_traits: list[str]
    extra_trait_taxa: list[str]
    pruned_missing_value_taxa: list[str]
    warnings: list[str]


@dataclass(slots=True)
class DiscreteMkTransformProfileRow:
    """One likelihood-profile row for a transformed discrete Mk fit."""

    transform_parameter_value: float
    log_likelihood: float


@dataclass(slots=True)
class DiscreteMkTransformWarning:
    """One transform-identifiability warning for a discrete Mk fit."""

    kind: str
    message: str


@dataclass(slots=True)
class DiscreteMkTransformFit:
    """One fitted transform surface layered over a discrete Mk likelihood."""

    transform_name: str
    parameter_name: str
    parameter_value: float
    lower_bound: float
    upper_bound: float
    starting_parameter_policy: str
    starting_parameter_value: float
    starting_parameter_log_likelihood: float
    coarse_grid_point_count: int
    fine_grid_point_count: int
    refinement_start_count: int
    function_evaluation_count: int
    hit_lower_parameter_boundary: bool
    hit_upper_parameter_boundary: bool
    transformed_tree_is_ultrametric: bool
    transformed_tree_minimum_tip_depth: float
    transformed_tree_maximum_tip_depth: float
    profile_rows: list[DiscreteMkTransformProfileRow]
    warnings: list[DiscreteMkTransformWarning]


@dataclass(slots=True)
class DiscreteMkTransformBaselineComparison:
    """Likelihood comparison against the untransformed branch-length surface."""

    baseline_transform: str
    baseline_log_likelihood: float
    baseline_parameter_count: int
    baseline_aic: float
    delta_log_likelihood: float
    delta_aic: float
    preferred_transform_by_aic: str


@dataclass(slots=True)
class DiscreteMkPatternLikelihoodRow:
    """One discrete-trait pattern likelihood row for total-likelihood reconstruction."""

    pattern_id: str
    pattern_weight: int
    tip_states: tuple[str, ...]
    raw_log_likelihood: float
    ascertainment_conditioning_log_probability: float | None
    log_likelihood: float


@dataclass(slots=True)
class DiscreteMkFitReport:
    """Discrete Mk trait-evolution fit over one rooted tree."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    ascertainment_policy: str
    state_ordering: str
    state_order: list[str]
    taxon_count: int
    input_audit: DiscreteMkInputAudit
    log_likelihood: float
    ascertainment_conditioning_log_probability: float | None
    invariant_pattern_log_probability: float | None
    parameter_count: int
    aic: float
    aicc: float
    likelihood_constant_policy: str
    likelihood_comparison_policy: str
    pattern_likelihood_rows: list[DiscreteMkPatternLikelihoodRow]
    transition_rate_rows: list[DiscreteTransitionRateRow]
    allowed_transition_pairs: list[tuple[str, str]]
    optimizer_diagnostics: DiscreteOptimizerDiagnostics
    overparameterized: bool
    transform_fit: DiscreteMkTransformFit | None
    transform_baseline_comparison: DiscreteMkTransformBaselineComparison | None
    baseline_comparison: DiscreteModelBaselineComparison | None


@dataclass(slots=True)
class DiscreteMkModelComparisonReport:
    """AIC/AICc model-comparison surface over governed discrete Mk fits."""

    tree_path: Path
    traits_path: Path
    trait: str
    ascertainment_policy: str
    taxon_count: int
    rows: list[ComparativeModelComparisonRow]
    better_model: str
    likelihood_constant_policy: str | None
    likelihood_comparison_policy: str
    noncomparable_likelihood_models: list[str]
    model_confidence_weight_basis: str
    model_confidence_delta_threshold: float
    selected_model_akaike_weight: float | None
    models_within_delta_aic_threshold: list[str]
    models_within_delta_aicc_threshold: list[str]
    uncertainty_language: str
    warnings: list[str]
