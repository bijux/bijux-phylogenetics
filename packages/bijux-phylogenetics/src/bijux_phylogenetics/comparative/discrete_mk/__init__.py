from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    AncestralDiscreteDataset,
    load_discrete_dataset,
    write_ancestral_rows,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteModelBaselineComparison,
    DiscreteOptimizerDiagnostics,
    DiscreteTransitionRateRow,
)
from bijux_phylogenetics.ancestral.discrete.likelihood import (
    build_transition_rate_rows as _build_transition_rate_rows,
    fit_discrete_mk_model as _fit_discrete_mk_model,
    tree_log_likelihood as _tree_log_likelihood,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    parameter_count as _parameter_count,
    resolve_discrete_model_name as _resolve_discrete_model_name,
    resolve_allowed_transition_pairs as _resolve_allowed_transition_pairs,
    resolve_state_order as _resolve_state_order,
)
from bijux_phylogenetics.ancestral.discrete.reconstruction import (
    _detect_discrete_overparameterization,
)
from bijux_phylogenetics.comparative.bounded_search import (
    BoundedSearchControls,
    run_bounded_maximization,
)
from bijux_phylogenetics.comparative.common import tip_root_depths
from bijux_phylogenetics.comparative.information_criteria import (
    compute_aic,
    compute_aicc,
    rank_model_comparison_rows,
)
from bijux_phylogenetics.comparative.evolutionary_modes import (
    transform_tree_for_evolutionary_mode,
)
from bijux_phylogenetics.comparative.models import ComparativeModelComparisonRow
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError
from bijux_phylogenetics.core.ultrametric import summarize_ultrametric_tip_depths

_DISCRETE_TRANSFORM_COARSE_GRID_POINT_COUNT = 9
_DISCRETE_TRANSFORM_FINE_GRID_POINT_COUNT = 17
_DISCRETE_TRANSFORM_REFINEMENT_START_COUNT = 1
_DISCRETE_DELTA_LOWER_BOUND = math.exp(-5.0)
_DISCRETE_DELTA_UPPER_BOUND = 3.0
_DISCRETE_EARLY_BURST_LOWER_BOUND = -10.0
_DISCRETE_EARLY_BURST_UPPER_BOUND = 10.0
DISCRETE_MK_MODEL_COMPARISON_ORDER = (
    "equal-rates",
    "symmetric",
    "all-rates-different",
)
DISCRETE_MK_LIKELIHOOD_CONSTANT_POLICY = (
    "continuous-time-markov-pruning-loglikelihood-has-no-extra-normalizing-constant"
)
DISCRETE_MK_LIKELIHOOD_COMPARISON_POLICY = (
    "raw-loglikelihood-and-derived-aic-are-directly-comparable-when-all-candidate-mk-models-share-the-owned-pruning-likelihood-policy"
)
DISCRETE_MK_MODEL_RANKING_POLICY = (
    "relative-aic-and-aicc-ranking-is-permitted-only-when-all-candidate-discrete-mk-models-share-one-pruning-likelihood-policy"
)
DISCRETE_MK_MODEL_CONFIDENCE_WEIGHT_BASIS = "AICc"
DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD = 2.0


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
class DiscreteMkFitReport:
    """Discrete Mk trait-evolution fit over one rooted tree."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    state_ordering: str
    state_order: list[str]
    taxon_count: int
    input_audit: DiscreteMkInputAudit
    log_likelihood: float
    parameter_count: int
    aic: float
    aicc: float
    likelihood_constant_policy: str
    likelihood_comparison_policy: str
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


def _normalize_transition_rate_rows(
    rows: list[DiscreteTransitionRateRow],
    *,
    state_ordering: str,
) -> list[DiscreteTransitionRateRow]:
    if state_ordering != "unordered":
        return rows
    return [
        DiscreteTransitionRateRow(
            source_state=row.source_state,
            target_state=row.target_state,
            transition_allowed=row.transition_allowed,
            step_distance=(1 if row.transition_allowed else row.step_distance),
            rate=row.rate,
        )
        for row in rows
    ]


def fit_discrete_mk_model(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    transform: str | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 1.0),
    delta_bounds: tuple[float, float] = (
        _DISCRETE_DELTA_LOWER_BOUND,
        _DISCRETE_DELTA_UPPER_BOUND,
    ),
    early_burst_bounds: tuple[float, float] = (
        _DISCRETE_EARLY_BURST_LOWER_BOUND,
        _DISCRETE_EARLY_BURST_UPPER_BOUND,
    ),
) -> DiscreteMkFitReport:
    """Fit one Mk discrete-trait model on a rooted tree."""
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    return fit_discrete_mk_model_from_dataset(
        dataset,
        model=model,
        transform=transform,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        allowed_transition_pairs=allowed_transition_pairs,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        early_burst_bounds=early_burst_bounds,
    )


def fit_discrete_mk_model_from_dataset(
    dataset: AncestralDiscreteDataset,
    *,
    model: str = "equal-rates",
    transform: str | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 1.0),
    delta_bounds: tuple[float, float] = (
        _DISCRETE_DELTA_LOWER_BOUND,
        _DISCRETE_DELTA_UPPER_BOUND,
    ),
    early_burst_bounds: tuple[float, float] = (
        _DISCRETE_EARLY_BURST_LOWER_BOUND,
        _DISCRETE_EARLY_BURST_UPPER_BOUND,
    ),
) -> DiscreteMkFitReport:
    """Fit one Mk discrete-trait model from a native discrete dataset."""
    resolved_model = _resolve_discrete_model_name(model)
    resolved_transform = _resolve_discrete_transform_name(transform)
    state_order = _resolve_state_order(
        dataset.observed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    _validate_discrete_transform_request(
        dataset,
        transform=resolved_transform,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        early_burst_bounds=early_burst_bounds,
    )
    resolved_allowed_transition_pairs = _resolve_allowed_transition_pairs(
        state_order,
        model=resolved_model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    (
        fit_tree,
        rate_matrix,
        root_prior,
        optimizer_diagnostics,
        transform_fit,
        transform_warning_rows,
    ) = _fit_discrete_mk_surface(
        dataset,
        model=resolved_model,
        transform=resolved_transform,
        state_ordering=state_ordering,
        state_order=state_order,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        early_burst_bounds=early_burst_bounds,
    )
    log_likelihood = _tree_log_likelihood(
        fit_tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        root_prior_mode="observed",
    )
    parameter_count = _discrete_parameter_count(
        state_count=len(state_order),
        model=resolved_model,
        transform=resolved_transform,
        state_ordering=state_ordering,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
    )
    aic = compute_aic(log_likelihood, parameter_count=parameter_count)
    aicc = compute_aicc(
        aic,
        sample_size=len(dataset.taxa),
        parameter_count=parameter_count,
    )
    baseline_comparison: DiscreteModelBaselineComparison | None = None
    transform_baseline_comparison: DiscreteMkTransformBaselineComparison | None = None
    if (
        resolved_model != "equal-rates"
        and allowed_transition_pairs is None
        and state_ordering == "unordered"
    ):
        baseline_fit = fit_discrete_mk_model_from_dataset(
            dataset,
            model="equal-rates",
            transform=resolved_transform,
            state_ordering=state_ordering,
            ordered_states=state_order,
            allowed_transition_pairs=None,
            lambda_bounds=lambda_bounds,
            kappa_bounds=kappa_bounds,
            delta_bounds=delta_bounds,
            early_burst_bounds=early_burst_bounds,
        )
        baseline_comparison = DiscreteModelBaselineComparison(
            baseline_model="equal-rates",
            baseline_log_likelihood=baseline_fit.log_likelihood,
            baseline_parameter_count=baseline_fit.parameter_count,
            baseline_aic=baseline_fit.aic,
            delta_log_likelihood=log_likelihood - baseline_fit.log_likelihood,
            delta_aic=aic - baseline_fit.aic,
            preferred_model_by_aic=(
                resolved_model if aic <= baseline_fit.aic else "equal-rates"
            ),
        )
    if resolved_transform is not None:
        transform_baseline_fit = fit_discrete_mk_model_from_dataset(
            dataset,
            model=resolved_model,
            transform=None,
            state_ordering=state_ordering,
            ordered_states=state_order,
            allowed_transition_pairs=allowed_transition_pairs,
            lambda_bounds=lambda_bounds,
            kappa_bounds=kappa_bounds,
            delta_bounds=delta_bounds,
            early_burst_bounds=early_burst_bounds,
        )
        transform_baseline_comparison = DiscreteMkTransformBaselineComparison(
            baseline_transform="untransformed",
            baseline_log_likelihood=transform_baseline_fit.log_likelihood,
            baseline_parameter_count=transform_baseline_fit.parameter_count,
            baseline_aic=transform_baseline_fit.aic,
            delta_log_likelihood=log_likelihood - transform_baseline_fit.log_likelihood,
            delta_aic=aic - transform_baseline_fit.aic,
            preferred_transform_by_aic=(
                resolved_transform if aic <= transform_baseline_fit.aic else "untransformed"
            ),
        )
    overparameterized = _detect_discrete_overparameterization(
        taxon_count=len(dataset.taxa),
        parameter_count=parameter_count,
    )
    warnings = list(dataset.warnings)
    if overparameterized:
        warnings.append(
            "the discrete Mk likelihood fit is likely overparameterized relative to the analyzed taxon count"
        )
    for warning in transform_warning_rows:
        warnings.append(warning.message)
    if not optimizer_diagnostics.converged:
        warnings.append(
            "the discrete Mk optimizer did not converge and should be interpreted cautiously"
        )
    if (
        optimizer_diagnostics.hit_lower_parameter_bound
        or optimizer_diagnostics.hit_upper_parameter_bound
    ):
        warnings.append(
            "one or more discrete Mk rate parameters hit an optimizer bound and should be interpreted as weakly identified"
        )
    if (
        baseline_comparison is not None
        and baseline_comparison.preferred_model_by_aic == "equal-rates"
    ):
        warnings.append(
            "the equal-rates baseline remains preferred by AIC over the requested discrete Mk model"
        )
    if (
        transform_baseline_comparison is not None
        and transform_baseline_comparison.preferred_transform_by_aic == "untransformed"
    ):
        warnings.append(
            "the untransformed branch-length baseline remains preferred by AIC over the requested discrete Mk transform"
        )
    input_audit = _build_discrete_mk_input_audit(dataset, warnings=warnings)
    return DiscreteMkFitReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        model=resolved_model,
        state_ordering=state_ordering,
        state_order=state_order,
        taxon_count=len(dataset.taxa),
        input_audit=input_audit,
        log_likelihood=log_likelihood,
        parameter_count=parameter_count,
        aic=aic,
        aicc=aicc,
        likelihood_constant_policy=DISCRETE_MK_LIKELIHOOD_CONSTANT_POLICY,
        likelihood_comparison_policy=DISCRETE_MK_LIKELIHOOD_COMPARISON_POLICY,
        transition_rate_rows=_normalize_transition_rate_rows(
            _build_transition_rate_rows(
                state_order=state_order,
                state_ordering=state_ordering,
                rate_matrix=rate_matrix,
                allowed_transition_pairs=resolved_allowed_transition_pairs,
            ),
            state_ordering=state_ordering,
        ),
        allowed_transition_pairs=[
            (state_order[left_index], state_order[right_index])
            for left_index, right_index in sorted(resolved_allowed_transition_pairs)
        ],
        optimizer_diagnostics=optimizer_diagnostics,
        overparameterized=overparameterized,
        transform_fit=transform_fit,
        transform_baseline_comparison=transform_baseline_comparison,
        baseline_comparison=baseline_comparison,
    )


def write_discrete_mk_summary_table(path: Path, report: DiscreteMkFitReport) -> Path:
    """Write one flat summary ledger for a discrete Mk fit."""
    baseline = report.baseline_comparison
    diagnostics = report.optimizer_diagnostics
    transform_fit = report.transform_fit
    transform_baseline = report.transform_baseline_comparison
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "transform",
            "transform_parameter_name",
            "transform_parameter_value",
            "transform_lower_bound",
            "transform_upper_bound",
            "transform_starting_parameter_policy",
            "transform_starting_parameter_value",
            "transform_starting_parameter_log_likelihood",
            "transform_coarse_grid_point_count",
            "transform_fine_grid_point_count",
            "transform_refinement_start_count",
            "transform_function_evaluation_count",
            "transform_hit_lower_parameter_boundary",
            "transform_hit_upper_parameter_boundary",
            "transform_warning_count",
            "transform_tree_is_ultrametric",
            "transform_tree_minimum_tip_depth",
            "transform_tree_maximum_tip_depth",
            "transform_baseline",
            "transform_baseline_parameter_count",
            "transform_baseline_log_likelihood",
            "transform_baseline_aic",
            "transform_delta_log_likelihood",
            "transform_delta_aic",
            "preferred_transform_by_aic",
            "state_ordering",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_state_count",
            "sparse_state_count",
            "log_likelihood",
            "parameter_count",
            "aic",
            "aicc",
            "likelihood_constant_policy",
            "likelihood_comparison_policy",
            "optimizer_name",
            "optimizer_converged",
            "optimizer_iteration_count",
            "optimizer_function_evaluation_count",
            "optimizer_simplex_shrink_count",
            "optimizer_initial_candidate_count",
            "optimizer_best_initial_scale",
            "optimizer_hit_lower_parameter_bound",
            "optimizer_hit_upper_parameter_bound",
            "overparameterized",
            "baseline_model",
            "baseline_parameter_count",
            "baseline_log_likelihood",
            "baseline_aic",
            "delta_log_likelihood",
            "delta_aic",
            "preferred_model_by_aic",
        ],
        rows=[
            {
                "trait": report.trait,
                "taxon_column": report.taxon_column,
                "model": report.model,
                "transform": "" if transform_fit is None else transform_fit.transform_name,
                "transform_parameter_name": (
                    "" if transform_fit is None else transform_fit.parameter_name
                ),
                "transform_parameter_value": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.parameter_value, ".15g")
                ),
                "transform_lower_bound": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.lower_bound, ".15g")
                ),
                "transform_upper_bound": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.upper_bound, ".15g")
                ),
                "transform_starting_parameter_policy": (
                    ""
                    if transform_fit is None
                    else transform_fit.starting_parameter_policy
                ),
                "transform_starting_parameter_value": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.starting_parameter_value, ".15g")
                ),
                "transform_starting_parameter_log_likelihood": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.starting_parameter_log_likelihood, ".15g")
                ),
                "transform_coarse_grid_point_count": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.coarse_grid_point_count)
                ),
                "transform_fine_grid_point_count": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.fine_grid_point_count)
                ),
                "transform_refinement_start_count": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.refinement_start_count)
                ),
                "transform_function_evaluation_count": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.function_evaluation_count)
                ),
                "transform_hit_lower_parameter_boundary": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.hit_lower_parameter_boundary).lower()
                ),
                "transform_hit_upper_parameter_boundary": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.hit_upper_parameter_boundary).lower()
                ),
                "transform_warning_count": (
                    ""
                    if transform_fit is None
                    else str(len(transform_fit.warnings))
                ),
                "transform_tree_is_ultrametric": (
                    ""
                    if transform_fit is None
                    else str(transform_fit.transformed_tree_is_ultrametric).lower()
                ),
                "transform_tree_minimum_tip_depth": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.transformed_tree_minimum_tip_depth, ".15g")
                ),
                "transform_tree_maximum_tip_depth": (
                    ""
                    if transform_fit is None
                    else format(transform_fit.transformed_tree_maximum_tip_depth, ".15g")
                ),
                "transform_baseline": (
                    ""
                    if transform_baseline is None
                    else transform_baseline.baseline_transform
                ),
                "transform_baseline_parameter_count": ""
                if transform_baseline is None
                else str(transform_baseline.baseline_parameter_count),
                "transform_baseline_log_likelihood": ""
                if transform_baseline is None
                else format(transform_baseline.baseline_log_likelihood, ".15g"),
                "transform_baseline_aic": ""
                if transform_baseline is None
                else format(transform_baseline.baseline_aic, ".15g"),
                "transform_delta_log_likelihood": ""
                if transform_baseline is None
                else format(transform_baseline.delta_log_likelihood, ".15g"),
                "transform_delta_aic": ""
                if transform_baseline is None
                else format(transform_baseline.delta_aic, ".15g"),
                "preferred_transform_by_aic": ""
                if transform_baseline is None
                else transform_baseline.preferred_transform_by_aic,
                "state_ordering": report.state_ordering,
                "analyzed_taxon_count": str(report.taxon_count),
                "excluded_taxon_count": str(
                    len(report.input_audit.pruned_missing_value_taxa)
                ),
                "observed_state_count": str(len(report.input_audit.observed_states)),
                "sparse_state_count": str(len(report.input_audit.sparse_states)),
                "log_likelihood": format(report.log_likelihood, ".15g"),
                "parameter_count": str(report.parameter_count),
                "aic": format(report.aic, ".15g"),
                "aicc": "inf"
                if math.isinf(report.aicc)
                else format(report.aicc, ".15g"),
                "likelihood_constant_policy": report.likelihood_constant_policy,
                "likelihood_comparison_policy": report.likelihood_comparison_policy,
                "optimizer_name": diagnostics.optimizer_name,
                "optimizer_converged": str(diagnostics.converged).lower(),
                "optimizer_iteration_count": str(diagnostics.iteration_count),
                "optimizer_function_evaluation_count": str(
                    diagnostics.function_evaluation_count
                ),
                "optimizer_simplex_shrink_count": str(diagnostics.simplex_shrink_count),
                "optimizer_initial_candidate_count": str(
                    diagnostics.initial_candidate_count
                ),
                "optimizer_best_initial_scale": format(
                    diagnostics.best_initial_scale,
                    ".15g",
                ),
                "optimizer_hit_lower_parameter_bound": str(
                    diagnostics.hit_lower_parameter_bound
                ).lower(),
                "optimizer_hit_upper_parameter_bound": str(
                    diagnostics.hit_upper_parameter_bound
                ).lower(),
                "overparameterized": str(report.overparameterized).lower(),
                "baseline_model": "" if baseline is None else baseline.baseline_model,
                "baseline_parameter_count": ""
                if baseline is None
                else str(baseline.baseline_parameter_count),
                "baseline_log_likelihood": ""
                if baseline is None
                else format(baseline.baseline_log_likelihood, ".15g"),
                "baseline_aic": ""
                if baseline is None
                else format(baseline.baseline_aic, ".15g"),
                "delta_log_likelihood": ""
                if baseline is None
                else format(baseline.delta_log_likelihood, ".15g"),
                "delta_aic": ""
                if baseline is None
                else format(baseline.delta_aic, ".15g"),
                "preferred_model_by_aic": (
                    "" if baseline is None else baseline.preferred_model_by_aic
                ),
            }
        ],
    )


def write_discrete_mk_rate_table(path: Path, report: DiscreteMkFitReport) -> Path:
    """Write one directed rate-matrix ledger for a discrete Mk fit."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_state",
            "target_state",
            "transition_allowed",
            "step_distance",
            "rate",
        ],
        rows=[
            {
                "source_state": row.source_state,
                "target_state": row.target_state,
                "transition_allowed": str(row.transition_allowed).lower(),
                "step_distance": str(row.step_distance),
                "rate": format(row.rate, ".15g"),
            }
            for row in report.transition_rate_rows
        ],
    )


def compare_discrete_mk_model_ranking(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    models: tuple[str, ...] | None = None,
    transform: str | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 1.0),
    delta_bounds: tuple[float, float] = (
        _DISCRETE_DELTA_LOWER_BOUND,
        _DISCRETE_DELTA_UPPER_BOUND,
    ),
    early_burst_bounds: tuple[float, float] = (
        _DISCRETE_EARLY_BURST_LOWER_BOUND,
        _DISCRETE_EARLY_BURST_UPPER_BOUND,
    ),
) -> DiscreteMkModelComparisonReport:
    """Compare governed discrete Mk models by AIC and AICc."""
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    return compare_discrete_mk_model_ranking_from_dataset(
        dataset,
        models=models,
        transform=transform,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        allowed_transition_pairs=allowed_transition_pairs,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        early_burst_bounds=early_burst_bounds,
    )


def compare_discrete_mk_model_ranking_from_dataset(
    dataset: AncestralDiscreteDataset,
    *,
    models: tuple[str, ...] | None = None,
    transform: str | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 1.0),
    delta_bounds: tuple[float, float] = (
        _DISCRETE_DELTA_LOWER_BOUND,
        _DISCRETE_DELTA_UPPER_BOUND,
    ),
    early_burst_bounds: tuple[float, float] = (
        _DISCRETE_EARLY_BURST_LOWER_BOUND,
        _DISCRETE_EARLY_BURST_UPPER_BOUND,
    ),
) -> DiscreteMkModelComparisonReport:
    """Compare a selected governed discrete Mk model set by information criteria."""
    selected_models = _comparison_models(models)
    state_order = _resolve_state_order(
        dataset.observed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    rows: list[ComparativeModelComparisonRow] = []
    fits: dict[str, DiscreteMkFitReport] = {}
    comparison_warnings: list[str] = []
    for model in selected_models:
        resolved_model = _resolve_discrete_model_name(model)
        try:
            fit = fit_discrete_mk_model_from_dataset(
                dataset,
                model=resolved_model,
                transform=transform,
                state_ordering=state_ordering,
                ordered_states=state_order,
                allowed_transition_pairs=allowed_transition_pairs,
                lambda_bounds=lambda_bounds,
                kappa_bounds=kappa_bounds,
                delta_bounds=delta_bounds,
                early_burst_bounds=early_burst_bounds,
            )
        except ComparativeMethodError as error:
            rows.append(
                ComparativeModelComparisonRow(
                    model=resolved_model,
                    parameter_count=_comparison_parameter_count(
                        state_order=state_order,
                        model=resolved_model,
                        transform=transform,
                        state_ordering=state_ordering,
                        allowed_transition_pairs=allowed_transition_pairs,
                    ),
                    log_likelihood=math.nan,
                    aic=math.inf,
                    aicc=math.inf,
                    comparable=False,
                    comparability_note=str(error),
                    selected=False,
                    likelihood_constant_policy=DISCRETE_MK_LIKELIHOOD_CONSTANT_POLICY,
                )
            )
            comparison_warnings.append(
                f"{resolved_model} is not comparable on this dataset because the owned fit failed: {error}"
            )
            continue
        fits[resolved_model] = fit
        row = ComparativeModelComparisonRow(
            model=fit.model,
            parameter_count=fit.parameter_count,
            log_likelihood=fit.log_likelihood,
            aic=fit.aic,
            aicc=fit.aicc,
            likelihood_constant_policy=fit.likelihood_constant_policy,
        )
        if not math.isfinite(row.aicc):
            row.comparable = False
            row.comparability_note = (
                "sample size is too small to compute finite AICc for this parameter count"
            )
            comparison_warnings.append(
                f"{resolved_model} is not comparable on AICc because the retained taxon count is too small for a {row.parameter_count}-parameter fit"
            )
        rows.append(row)
    if not fits:
        raise ComparativeMethodError(
            "no discrete Mk model remained comparable for the requested dataset"
        )
    likelihood_constant_policy, noncomparable_likelihood_models = (
        rank_model_comparison_rows(
            rows,
            delta_threshold=DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        )
    )
    if noncomparable_likelihood_models:
        blocked_models = ", ".join(noncomparable_likelihood_models)
        comparison_warnings.append(
            "discrete Mk ranking excluded models with incompatible likelihood "
            f"constant policies: {blocked_models}"
        )
    selected_rows = [row for row in rows if row.selected]
    if not selected_rows:
        if noncomparable_likelihood_models:
            raise ComparativeMethodError(
                "mixed likelihood constant policies prevent ranking incompatible discrete Mk models"
            )
        raise ComparativeMethodError(
            "no finite AICc model remained available for discrete Mk comparison"
        )
    better_model = selected_rows[0].model
    if len(selected_rows) > 1:
        tied_models = ", ".join(row.model for row in selected_rows)
        comparison_warnings.append(
            f"multiple discrete Mk models remain tied at the selected AICc boundary: {tied_models}"
        )
    selected_fit = fits[better_model]
    if selected_fit.overparameterized:
        comparison_warnings.append(
            "selected discrete Mk fit remains overparameterized relative to the analyzed taxon count and should be interpreted cautiously"
        )
    if not selected_fit.optimizer_diagnostics.converged:
        comparison_warnings.append(
            "selected discrete Mk fit did not converge cleanly, so model-ranking confidence is reduced"
        )
    if (
        selected_fit.optimizer_diagnostics.hit_lower_parameter_bound
        or selected_fit.optimizer_diagnostics.hit_upper_parameter_bound
    ):
        comparison_warnings.append(
            "selected discrete Mk fit hits one or more optimizer bounds, so the winning rate surface should be treated as weakly identified"
        )
    return DiscreteMkModelComparisonReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        trait=dataset.trait,
        taxon_count=len(dataset.taxa),
        rows=rows,
        better_model=better_model,
        likelihood_constant_policy=likelihood_constant_policy,
        likelihood_comparison_policy=DISCRETE_MK_MODEL_RANKING_POLICY,
        noncomparable_likelihood_models=noncomparable_likelihood_models,
        model_confidence_weight_basis=DISCRETE_MK_MODEL_CONFIDENCE_WEIGHT_BASIS,
        model_confidence_delta_threshold=DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        selected_model_akaike_weight=_selected_model_akaike_weight(rows),
        models_within_delta_aic_threshold=_models_within_delta_threshold(
            rows,
            criterion="aic",
            threshold=DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        models_within_delta_aicc_threshold=_models_within_delta_threshold(
            rows,
            criterion="aicc",
            threshold=DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        uncertainty_language=_model_confidence_uncertainty_language(
            rows,
            better_model=better_model,
            threshold=DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        warnings=comparison_warnings,
    )


def _resolve_discrete_transform_name(transform: str | None) -> str | None:
    if transform is None:
        return None
    aliases = {
        "lambda": "lambda",
        "pagel-lambda": "lambda",
        "kappa": "kappa",
        "pagel-kappa": "kappa",
        "delta": "delta",
        "pagel-delta": "delta",
        "EB": "early-burst",
        "early-burst": "early-burst",
    }
    resolved = aliases.get(transform)
    if resolved is None:
        raise ComparativeMethodError(
            "unsupported discrete Mk transform; expected one of: lambda, kappa, delta, early-burst"
        )
    return resolved


def _validate_discrete_transform_request(
    dataset: AncestralDiscreteDataset,
    *,
    transform: str | None,
    state_ordering: str,
    allowed_transition_pairs: list[tuple[str, str]] | None,
    lambda_bounds: tuple[float, float],
    kappa_bounds: tuple[float, float],
    delta_bounds: tuple[float, float],
    early_burst_bounds: tuple[float, float],
) -> None:
    if transform is None:
        return
    if transform not in {"lambda", "kappa", "delta", "early-burst"}:
        raise ComparativeMethodError(
            "unsupported discrete Mk transform; expected one of: lambda, kappa, delta, early-burst"
        )
    if state_ordering != "unordered":
        raise ComparativeMethodError(
            "discrete Mk branch-length transforms currently support only unordered ER, SYM, and ARD fits"
        )
    if allowed_transition_pairs is not None:
        raise ComparativeMethodError(
            "discrete Mk branch-length transforms do not yet support custom transition constraints"
        )
    if transform == "lambda":
        ultrametric_summary = summarize_ultrametric_tip_depths(
            tip_root_depths(dataset.tree, dataset.taxa),
            tolerance=1e-12,
        )
        if not ultrametric_summary.ultrametric:
            raise ComparativeMethodError(
                "discrete Mk lambda transform requires an ultrametric rooted tree so the transformed branch lengths preserve tip-depth meaning"
        )
        _validate_lambda_bounds(lambda_bounds)
        return
    if transform == "kappa":
        _validate_kappa_bounds(kappa_bounds)
        return
    if transform == "delta":
        _validate_delta_bounds(delta_bounds)
        return
    _validate_early_burst_bounds(early_burst_bounds)


def _validate_lambda_bounds(bounds: tuple[float, float]) -> None:
    lower, upper = bounds
    if not 0.0 <= lower < upper <= 1.0:
        raise ComparativeMethodError(
            "discrete Mk lambda bounds must be strictly increasing within [0, 1]"
        )


def _validate_kappa_bounds(bounds: tuple[float, float]) -> None:
    lower, upper = bounds
    if not 0.0 <= lower < upper <= 1.0:
        raise ComparativeMethodError(
            "discrete Mk kappa bounds must be strictly increasing within [0, 1]"
        )


def _validate_delta_bounds(bounds: tuple[float, float]) -> None:
    lower, upper = bounds
    if not 0.0 < lower < upper <= _DISCRETE_DELTA_UPPER_BOUND:
        raise ComparativeMethodError(
            "discrete Mk delta bounds must be strictly increasing within "
            f"(0, {_DISCRETE_DELTA_UPPER_BOUND:g}]"
        )


def _validate_early_burst_bounds(bounds: tuple[float, float]) -> None:
    lower, upper = bounds
    if not (
        _DISCRETE_EARLY_BURST_LOWER_BOUND
        <= lower
        < upper
        <= _DISCRETE_EARLY_BURST_UPPER_BOUND
    ):
        raise ComparativeMethodError(
            "discrete Mk early-burst bounds must be strictly increasing within "
            f"[{_DISCRETE_EARLY_BURST_LOWER_BOUND:g}, {_DISCRETE_EARLY_BURST_UPPER_BOUND:g}]"
        )


def _comparison_models(models: tuple[str, ...] | None) -> tuple[str, ...]:
    if models is None:
        return DISCRETE_MK_MODEL_COMPARISON_ORDER
    deduplicated: list[str] = []
    unknown: list[str] = []
    for model in models:
        try:
            resolved_model = _resolve_discrete_model_name(model)
        except ValueError:
            unknown.append(model)
            continue
        if resolved_model not in DISCRETE_MK_MODEL_COMPARISON_ORDER:
            unknown.append(model)
            continue
        if resolved_model not in deduplicated:
            deduplicated.append(resolved_model)
    if unknown:
        raise ComparativeMethodError(
            "unsupported discrete Mk comparison model(s): "
            + ", ".join(sorted(set(unknown)))
        )
    if not deduplicated:
        raise ComparativeMethodError(
            "at least one supported discrete Mk model is required for comparison"
        )
    return tuple(deduplicated)


def _discrete_parameter_count(
    *,
    state_count: int,
    model: str,
    transform: str | None,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
) -> int:
    parameter_count = _parameter_count(
        state_count,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    if transform in {"lambda", "kappa", "delta", "early-burst"}:
        return parameter_count + 1
    return parameter_count


def _comparison_parameter_count(
    *,
    state_order: list[str],
    model: str,
    transform: str | None,
    state_ordering: str,
    allowed_transition_pairs: list[tuple[str, str]] | None,
) -> int:
    resolved_allowed_transition_pairs = _resolve_allowed_transition_pairs(
        state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    return _discrete_parameter_count(
        state_count=len(state_order),
        model=model,
        transform=_resolve_discrete_transform_name(transform),
        state_ordering=state_ordering,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
    )


def _fit_discrete_mk_surface(
    dataset: AncestralDiscreteDataset,
    *,
    model: str,
    transform: str | None,
    state_ordering: str,
    state_order: list[str],
    allowed_transition_pairs: set[tuple[int, int]],
    lambda_bounds: tuple[float, float],
    kappa_bounds: tuple[float, float],
    delta_bounds: tuple[float, float],
    early_burst_bounds: tuple[float, float],
) -> tuple[
    PhyloTree,
    object,
    object,
    DiscreteOptimizerDiagnostics,
    DiscreteMkTransformFit | None,
    list[DiscreteMkTransformWarning],
]:
    if transform is None:
        rate_matrix, root_prior, optimizer_diagnostics = _fit_discrete_mk_model(
            dataset.tree,
            dataset.states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            root_prior_mode="observed",
        )
        return (
            dataset.tree,
            rate_matrix,
            root_prior,
            optimizer_diagnostics,
            None,
            [],
        )
    if transform not in {"lambda", "kappa", "delta", "early-burst"}:
        raise ComparativeMethodError(
            "unsupported discrete Mk transform; expected one of: lambda, kappa, delta, early-burst"
        )
    if transform == "delta":
        bounds = delta_bounds
    elif transform == "early-burst":
        bounds = early_burst_bounds
    elif transform == "lambda":
        bounds = lambda_bounds
    else:
        bounds = kappa_bounds
    return _fit_discrete_mk_parameterized_transform_surface(
        dataset,
        model=model,
        transform=transform,
        state_ordering=state_ordering,
        state_order=state_order,
        allowed_transition_pairs=allowed_transition_pairs,
        bounds=bounds,
    )


def _fit_discrete_mk_parameterized_transform_surface(
    dataset: AncestralDiscreteDataset,
    *,
    model: str,
    transform: str,
    state_ordering: str,
    state_order: list[str],
    allowed_transition_pairs: set[tuple[int, int]],
    bounds: tuple[float, float],
) -> tuple[
    PhyloTree,
    object,
    object,
    DiscreteOptimizerDiagnostics,
    DiscreteMkTransformFit,
    list[DiscreteMkTransformWarning],
]:
    lower, upper = bounds
    search_result = run_bounded_maximization(
        lower_bound=lower,
        upper_bound=upper,
        controls=BoundedSearchControls(
            coarse_grid_point_count=_DISCRETE_TRANSFORM_COARSE_GRID_POINT_COUNT,
            fine_grid_point_count=_DISCRETE_TRANSFORM_FINE_GRID_POINT_COUNT,
            refinement_start_count=_DISCRETE_TRANSFORM_REFINEMENT_START_COUNT,
        ),
        evaluate=lambda parameter_value: _evaluate_discrete_mk_transform_search_candidate(
            dataset,
            model=model,
            transform=transform,
            state_ordering=state_ordering,
            state_order=state_order,
            allowed_transition_pairs=allowed_transition_pairs,
            parameter_value=parameter_value,
        ),
        optimizer_name="governed-two-stage-grid-search",
        parameter_search_strategy="bounded-two-stage-grid-search",
    )
    transformed_tree, rate_matrix, root_prior, optimizer_diagnostics = (
        search_result.payload
    )
    best_parameter_value = search_result.parameter_value
    transform_warning_rows = _discrete_transform_warning_rows(
        transform=transform,
        parameter_value=best_parameter_value,
        profile=[
            (row.parameter_value, row.objective_value)
            for row in search_result.profile_rows
        ],
        bounds=bounds,
    )
    transformed_ultrametric = summarize_ultrametric_tip_depths(
        tip_root_depths(transformed_tree, dataset.taxa),
        tolerance=1e-12,
    )
    transform_fit = DiscreteMkTransformFit(
        transform_name=transform,
        parameter_name=_discrete_transform_parameter_name(transform),
        parameter_value=best_parameter_value,
        lower_bound=lower,
        upper_bound=upper,
        starting_parameter_policy=search_result.diagnostics.starting_parameter_policy,
        starting_parameter_value=search_result.diagnostics.starting_parameter_value,
        starting_parameter_log_likelihood=search_result.diagnostics.starting_parameter_objective_value,
        coarse_grid_point_count=search_result.diagnostics.coarse_grid_point_count,
        fine_grid_point_count=search_result.diagnostics.fine_grid_point_count,
        refinement_start_count=search_result.diagnostics.refinement_start_count,
        function_evaluation_count=search_result.diagnostics.function_evaluation_count,
        hit_lower_parameter_boundary=search_result.diagnostics.hit_lower_boundary,
        hit_upper_parameter_boundary=search_result.diagnostics.hit_upper_boundary,
        transformed_tree_is_ultrametric=transformed_ultrametric.ultrametric,
        transformed_tree_minimum_tip_depth=transformed_ultrametric.minimum_tip_depth,
        transformed_tree_maximum_tip_depth=transformed_ultrametric.maximum_tip_depth,
        profile_rows=[
            DiscreteMkTransformProfileRow(
                transform_parameter_value=row.parameter_value,
                log_likelihood=row.objective_value,
            )
            for row in search_result.profile_rows
        ],
        warnings=transform_warning_rows,
    )
    return (
        transformed_tree,
        rate_matrix,
        root_prior,
        optimizer_diagnostics,
        transform_fit,
        transform_warning_rows,
    )


def _evaluate_discrete_mk_transform_candidate(
    dataset: AncestralDiscreteDataset,
    *,
    model: str,
    transform: str,
    state_ordering: str,
    state_order: list[str],
    allowed_transition_pairs: set[tuple[int, int]],
    parameter_value: float,
) -> tuple[object, object, object, DiscreteOptimizerDiagnostics, float]:
    parameter_value = _transform_cache_key(parameter_value)
    transformed_tree = _transform_discrete_mk_tree(
        dataset.tree,
        transform=transform,
        parameter_value=parameter_value,
    )
    rate_matrix, root_prior, optimizer_diagnostics = _fit_discrete_mk_model(
        transformed_tree,
        dataset.states_by_taxon,
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
        root_prior_mode="observed",
    )
    log_likelihood = _tree_log_likelihood(
        transformed_tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        root_prior_mode="observed",
    )
    return (
        transformed_tree,
        rate_matrix,
        root_prior,
        optimizer_diagnostics,
        log_likelihood,
    )


def _evaluate_discrete_mk_transform_search_candidate(
    dataset: AncestralDiscreteDataset,
    *,
    model: str,
    transform: str,
    state_ordering: str,
    state_order: list[str],
    allowed_transition_pairs: set[tuple[int, int]],
    parameter_value: float,
) -> tuple[tuple[PhyloTree, object, object, DiscreteOptimizerDiagnostics], float]:
    (
        transformed_tree,
        rate_matrix,
        root_prior,
        optimizer_diagnostics,
        log_likelihood,
    ) = _evaluate_discrete_mk_transform_candidate(
        dataset,
        model=model,
        transform=transform,
        state_ordering=state_ordering,
        state_order=state_order,
        allowed_transition_pairs=allowed_transition_pairs,
        parameter_value=parameter_value,
    )
    return (
        transformed_tree,
        rate_matrix,
        root_prior,
        optimizer_diagnostics,
    ), log_likelihood


def _discrete_transform_mode_name(transform: str) -> str:
    if transform == "lambda":
        return "pagel-lambda"
    if transform == "kappa":
        return "pagel-kappa"
    if transform == "delta":
        return "pagel-delta"
    if transform == "early-burst":
        return "early-burst"
    raise ComparativeMethodError(
        "unsupported discrete Mk transform; expected one of: lambda, kappa, delta, early-burst"
    )


def _discrete_transform_parameter_name(transform: str) -> str:
    if transform == "early-burst":
        return "a"
    return transform


def _transform_discrete_mk_tree(
    tree: PhyloTree,
    *,
    transform: str,
    parameter_value: float,
) -> PhyloTree:
    transformed_parameter_value = (
        -parameter_value if transform == "early-burst" else parameter_value
    )
    return transform_tree_for_evolutionary_mode(
        tree,
        mode=_discrete_transform_mode_name(transform),
        parameter_value=transformed_parameter_value,
        sigsq=1.0,
    )


def _transform_cache_key(parameter_value: float) -> float:
    return float(format(parameter_value, ".15g"))


def _discrete_transform_warning_rows(
    *,
    transform: str,
    parameter_value: float,
    profile: list[tuple[float, float]],
    bounds: tuple[float, float],
) -> list[DiscreteMkTransformWarning]:
    lower, upper = bounds
    span = upper - lower
    boundary_tolerance = max(span / 160.0, 1e-9)
    ordered_log_likelihoods = sorted(
        (log_likelihood for _, log_likelihood in profile),
        reverse=True,
    )
    warnings: list[DiscreteMkTransformWarning] = []
    if math.isclose(parameter_value, lower, rel_tol=0.0, abs_tol=boundary_tolerance) or math.isclose(
        parameter_value, upper, rel_tol=0.0, abs_tol=boundary_tolerance
    ):
        warnings.append(
            DiscreteMkTransformWarning(
                kind=(
                    "boundary_early_burst"
                    if transform == "early-burst"
                    else f"boundary_{transform}"
                ),
                message=(
                    "best-supported discrete Mk early-burst transform parameter falls on the search boundary and may not be well identified"
                    if transform == "early-burst"
                    else f"best-supported discrete Mk {transform} falls on the search boundary and may not be well identified"
                ),
            )
        )
    if len(ordered_log_likelihoods) > 1 and (
        ordered_log_likelihoods[0] - ordered_log_likelihoods[1] < 0.5
    ):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="flat_likelihood",
                message=f"discrete Mk {transform} likelihood stays shallow across the bounded search, so branch-length transformation support may be unstable",
            )
        )
    if transform == "lambda" and parameter_value <= lower + max(boundary_tolerance, 1e-6):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="weak_phylogenetic_signal",
                message="best-supported discrete Mk lambda remains close to the zero-signal boundary and may be difficult to distinguish from a star-like transition surface",
            )
        )
    if transform == "lambda" and parameter_value >= upper - max(boundary_tolerance, 1e-6):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="brownian_limit",
                message="best-supported discrete Mk lambda remains close to the untransformed boundary and may be difficult to distinguish from the original branch-length surface",
            )
        )
    if transform == "kappa" and parameter_value <= lower + max(boundary_tolerance, 1e-6):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="branch_length_flattening_limit",
                message="best-supported discrete Mk kappa remains close to the zero-contrast boundary and may be difficult to distinguish from branch-length flattening",
            )
        )
    if transform == "kappa" and parameter_value >= upper - max(boundary_tolerance, 1e-6):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="branch_length_contrast_limit",
                message="best-supported discrete Mk kappa remains close to the untransformed branch-length boundary and may be weakly identified",
            )
        )
    if transform == "delta" and parameter_value <= lower + max(boundary_tolerance, 1e-6):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="earliest_change_limit",
                message="best-supported discrete Mk delta remains close to the earliest-change boundary and may be difficult to distinguish from an extreme root-concentrated transition surface",
            )
        )
    if transform == "delta" and parameter_value >= upper - max(boundary_tolerance, 1e-6):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="late_change_limit",
                message="best-supported discrete Mk delta remains close to the late-change boundary and may be difficult to distinguish from an extreme tip-concentrated transition surface",
            )
        )
    if transform == "early-burst" and abs(parameter_value) <= max(
        boundary_tolerance, 1e-6
    ):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="brownian_like_rate_change",
                message="best-supported discrete Mk early-burst transform parameter remains close to the zero-change boundary and may be difficult to distinguish from the untransformed branch-length surface",
            )
        )
    if transform == "early-burst" and parameter_value <= lower + max(
        boundary_tolerance, 1e-6
    ):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="late_change_limit",
                message="best-supported discrete Mk early-burst rate change remains close to the strongest late-change boundary and may be difficult to distinguish from an extreme tip-concentrated transition surface",
            )
        )
    if transform == "early-burst" and parameter_value >= upper - max(
        boundary_tolerance, 1e-6
    ):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="early_change_limit",
                message="best-supported discrete Mk early-burst rate change remains close to the strongest early-change boundary and may be difficult to distinguish from an extreme root-concentrated transition surface",
            )
        )
    return warnings


def _build_discrete_mk_input_audit(
    dataset: AncestralDiscreteDataset,
    *,
    warnings: list[str],
) -> DiscreteMkInputAudit:
    ultrametric_summary = summarize_ultrametric_tip_depths(
        tip_root_depths(dataset.tree, dataset.taxa),
        tolerance=1e-12,
    )
    return DiscreteMkInputAudit(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        trait=dataset.trait,
        taxon_count=len(dataset.taxa),
        taxa=list(dataset.taxa),
        observed_states=list(dataset.observed_states),
        state_counts=dict(dataset.state_counts),
        sparse_states=list(dataset.sparse_states),
        tree_is_ultrametric=ultrametric_summary.ultrametric,
        minimum_root_to_tip_depth=ultrametric_summary.minimum_tip_depth,
        maximum_root_to_tip_depth=ultrametric_summary.maximum_tip_depth,
        ultrametric_policy="accept-rooted-trees-and-report-ultrametricity",
        missing_value_policy="prune-overlapping-missing-values",
        missing_from_traits=list(dataset.alignment_report.dropped_tree_taxa),
        extra_trait_taxa=list(dataset.alignment_report.dropped_trait_taxa),
        pruned_missing_value_taxa=list(dataset.dropped_missing_taxa),
        warnings=warnings,
    )


def _selected_model_akaike_weight(
    rows: list[ComparativeModelComparisonRow],
) -> float | None:
    selected_row = next((row for row in rows if row.selected), None)
    if selected_row is None:
        return None
    return selected_row.akaike_weight


def _models_within_delta_threshold(
    rows: list[ComparativeModelComparisonRow],
    *,
    criterion: str,
    threshold: float,
) -> list[str]:
    selected_rows: list[str] = []
    for row in rows:
        if not row.comparable:
            continue
        within_threshold = (
            row.within_delta_aic_threshold
            if criterion == "aic"
            else row.within_delta_aicc_threshold
        )
        if within_threshold:
            selected_rows.append(row.model)
    return selected_rows


def _model_confidence_uncertainty_language(
    rows: list[ComparativeModelComparisonRow],
    *,
    better_model: str,
    threshold: float,
) -> str:
    selected_row = next((row for row in rows if row.model == better_model), None)
    if selected_row is None or selected_row.akaike_weight is None:
        return (
            "model confidence is unresolved because no comparable discrete Mk "
            "candidate retained one finite AICc surface for Akaike-weight review"
        )
    nearby_models = [
        row.model
        for row in rows
        if row.comparable
        and row.model != better_model
        and row.within_delta_aicc_threshold is True
    ]
    if nearby_models:
        tied_models = ", ".join([better_model, *nearby_models])
        return (
            "model confidence is limited because "
            f"{tied_models} remain within {threshold:.1f} AICc units of the selected "
            f"surface; {better_model} carries Akaike weight "
            f"{selected_row.akaike_weight:.3f}"
        )
    if selected_row.akaike_weight < 0.9:
        return (
            f"model confidence is moderate because {better_model} carries Akaike "
            f"weight {selected_row.akaike_weight:.3f} even though no runner-up remains "
            f"within {threshold:.1f} AICc units"
        )
    return (
        f"model confidence is strong because {better_model} carries Akaike weight "
        f"{selected_row.akaike_weight:.3f} and no runner-up remains within "
        f"{threshold:.1f} AICc units"
    )
