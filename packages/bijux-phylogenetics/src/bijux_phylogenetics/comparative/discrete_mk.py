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
    _build_transition_rate_rows,
    _detect_discrete_overparameterization,
    _fit_discrete_mk_model,
    _parameter_count,
    _resolve_allowed_transition_pairs,
    _resolve_discrete_model_name,
    _resolve_state_order,
    _tree_log_likelihood,
)
from bijux_phylogenetics.comparative.common import tip_root_depths
from bijux_phylogenetics.comparative.evolutionary_modes import (
    transform_tree_for_evolutionary_mode,
)
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError
from bijux_phylogenetics.core.ultrametric import summarize_ultrametric_tip_depths

_DISCRETE_TRANSFORM_COARSE_GRID_POINT_COUNT = 9
_DISCRETE_TRANSFORM_FINE_GRID_POINT_COUNT = 17


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
    coarse_grid_point_count: int
    fine_grid_point_count: int
    function_evaluation_count: int
    hit_lower_parameter_boundary: bool
    hit_upper_parameter_boundary: bool
    transformed_tree_is_ultrametric: bool
    transformed_tree_minimum_tip_depth: float
    transformed_tree_maximum_tip_depth: float
    profile_rows: list[DiscreteMkTransformProfileRow]
    warnings: list[DiscreteMkTransformWarning]


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
    transition_rate_rows: list[DiscreteTransitionRateRow]
    allowed_transition_pairs: list[tuple[str, str]]
    optimizer_diagnostics: DiscreteOptimizerDiagnostics
    overparameterized: bool
    transform_fit: DiscreteMkTransformFit | None
    baseline_comparison: DiscreteModelBaselineComparison | None


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
    aic = (2.0 * parameter_count) - (2.0 * log_likelihood)
    aicc = _fit_aicc(
        aic, sample_size=len(dataset.taxa), parameter_count=parameter_count
    )
    baseline_comparison: DiscreteModelBaselineComparison | None = None
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
        baseline_comparison=baseline_comparison,
    )


def write_discrete_mk_summary_table(path: Path, report: DiscreteMkFitReport) -> Path:
    """Write one flat summary ledger for a discrete Mk fit."""
    baseline = report.baseline_comparison
    diagnostics = report.optimizer_diagnostics
    transform_fit = report.transform_fit
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
            "transform_coarse_grid_point_count",
            "transform_fine_grid_point_count",
            "transform_function_evaluation_count",
            "transform_hit_lower_parameter_boundary",
            "transform_hit_upper_parameter_boundary",
            "transform_warning_count",
            "transform_tree_is_ultrametric",
            "transform_tree_minimum_tip_depth",
            "transform_tree_maximum_tip_depth",
            "state_ordering",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_state_count",
            "sparse_state_count",
            "log_likelihood",
            "parameter_count",
            "aic",
            "aicc",
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


def _resolve_discrete_transform_name(transform: str | None) -> str | None:
    if transform is None:
        return None
    aliases = {
        "lambda": "lambda",
        "pagel-lambda": "lambda",
        "kappa": "kappa",
        "pagel-kappa": "kappa",
    }
    resolved = aliases.get(transform)
    if resolved is None:
        raise ComparativeMethodError(
            "unsupported discrete Mk transform; expected one of: lambda, kappa"
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
) -> None:
    if transform is None:
        return
    if transform not in {"lambda", "kappa"}:
        raise ComparativeMethodError(
            "unsupported discrete Mk transform; expected one of: lambda, kappa"
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
    _validate_kappa_bounds(kappa_bounds)


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
    if transform in {"lambda", "kappa"}:
        return parameter_count + 1
    return parameter_count


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
    if transform not in {"lambda", "kappa"}:
        raise ComparativeMethodError(
            "unsupported discrete Mk transform; expected one of: lambda, kappa"
        )
    return _fit_discrete_mk_parameterized_transform_surface(
        dataset,
        model=model,
        transform=transform,
        state_ordering=state_ordering,
        state_order=state_order,
        allowed_transition_pairs=allowed_transition_pairs,
        bounds=lambda_bounds if transform == "lambda" else kappa_bounds,
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
    candidate_cache: dict[
        float, tuple[PhyloTree, object, object, DiscreteOptimizerDiagnostics, float]
    ] = {}

    def evaluate(
        parameter_value: float,
    ) -> tuple[PhyloTree, object, object, DiscreteOptimizerDiagnostics, float]:
        cache_key = _transform_cache_key(parameter_value)
        cached = candidate_cache.get(cache_key)
        if cached is not None:
            return cached
        candidate = _evaluate_discrete_mk_transform_candidate(
            dataset,
            model=model,
            transform=transform,
            state_ordering=state_ordering,
            state_order=state_order,
            allowed_transition_pairs=allowed_transition_pairs,
            parameter_value=cache_key,
        )
        candidate_cache[cache_key] = candidate
        return candidate

    coarse_candidates = _grid_values(
        lower,
        upper,
        point_count=_DISCRETE_TRANSFORM_COARSE_GRID_POINT_COUNT,
    )
    coarse_rows: list[tuple[float, float]] = []
    best_fit: tuple[
        float, PhyloTree, object, object, DiscreteOptimizerDiagnostics, float
    ] | None = None
    for parameter_value in coarse_candidates:
        (
            transformed_tree,
            rate_matrix,
            root_prior,
            optimizer_diagnostics,
            log_likelihood,
        ) = evaluate(parameter_value)
        coarse_rows.append((parameter_value, log_likelihood))
        if best_fit is None or log_likelihood > best_fit[-1]:
            best_fit = (
                parameter_value,
                transformed_tree,
                rate_matrix,
                root_prior,
                optimizer_diagnostics,
                log_likelihood,
            )
    assert best_fit is not None
    best_coarse_index = max(
        range(len(coarse_rows)),
        key=lambda index: coarse_rows[index][1],
    )
    fine_lower, fine_upper = _transform_search_bracket(
        coarse_candidates=coarse_candidates,
        best_index=best_coarse_index,
        lower=lower,
        upper=upper,
    )
    profile_rows: list[tuple[float, float]] = list(coarse_rows)
    if not math.isclose(fine_lower, fine_upper, rel_tol=0.0, abs_tol=1e-15):
        phi = (math.sqrt(5.0) - 1.0) / 2.0
        left = fine_upper - phi * (fine_upper - fine_lower)
        right = fine_lower + phi * (fine_upper - fine_lower)
        (
            left_tree,
            left_rate_matrix,
            left_root_prior,
            left_optimizer_diagnostics,
            left_score,
        ) = evaluate(left)
        profile_rows.append((left, left_score))
        if left_score > best_fit[-1]:
            best_fit = (
                left,
                left_tree,
                left_rate_matrix,
                left_root_prior,
                left_optimizer_diagnostics,
                left_score,
            )
        (
            right_tree,
            right_rate_matrix,
            right_root_prior,
            right_optimizer_diagnostics,
            right_score,
        ) = evaluate(right)
        profile_rows.append((right, right_score))
        if right_score > best_fit[-1]:
            best_fit = (
                right,
                right_tree,
                right_rate_matrix,
                right_root_prior,
                right_optimizer_diagnostics,
                right_score,
            )
        for _ in range(max(_DISCRETE_TRANSFORM_FINE_GRID_POINT_COUNT - 2, 0)):
            if abs(fine_upper - fine_lower) < 1e-6:
                break
            if left_score > right_score:
                fine_upper = right
                right = left
                right_tree = left_tree
                right_rate_matrix = left_rate_matrix
                right_root_prior = left_root_prior
                right_optimizer_diagnostics = left_optimizer_diagnostics
                right_score = left_score
                left = fine_upper - phi * (fine_upper - fine_lower)
                (
                    left_tree,
                    left_rate_matrix,
                    left_root_prior,
                    left_optimizer_diagnostics,
                    left_score,
                ) = evaluate(left)
                profile_rows.append((left, left_score))
                if left_score > best_fit[-1]:
                    best_fit = (
                        left,
                        left_tree,
                        left_rate_matrix,
                        left_root_prior,
                        left_optimizer_diagnostics,
                        left_score,
                    )
            else:
                fine_lower = left
                left = right
                left_tree = right_tree
                left_rate_matrix = right_rate_matrix
                left_root_prior = right_root_prior
                left_optimizer_diagnostics = right_optimizer_diagnostics
                left_score = right_score
                right = fine_lower + phi * (fine_upper - fine_lower)
                (
                    right_tree,
                    right_rate_matrix,
                    right_root_prior,
                    right_optimizer_diagnostics,
                    right_score,
                ) = evaluate(right)
                profile_rows.append((right, right_score))
                if right_score > best_fit[-1]:
                    best_fit = (
                        right,
                        right_tree,
                        right_rate_matrix,
                        right_root_prior,
                        right_optimizer_diagnostics,
                        right_score,
                    )
    (
        best_parameter_value,
        transformed_tree,
        rate_matrix,
        root_prior,
        optimizer_diagnostics,
        _,
    ) = best_fit
    transform_warning_rows = _discrete_transform_warning_rows(
        transform=transform,
        parameter_value=best_parameter_value,
        profile=profile_rows,
        bounds=bounds,
    )
    transformed_ultrametric = summarize_ultrametric_tip_depths(
        tip_root_depths(transformed_tree, dataset.taxa),
        tolerance=1e-12,
    )
    transform_fit = DiscreteMkTransformFit(
        transform_name=transform,
        parameter_name=transform,
        parameter_value=best_parameter_value,
        lower_bound=lower,
        upper_bound=upper,
        coarse_grid_point_count=_DISCRETE_TRANSFORM_COARSE_GRID_POINT_COUNT,
        fine_grid_point_count=_DISCRETE_TRANSFORM_FINE_GRID_POINT_COUNT,
        function_evaluation_count=len(candidate_cache),
        hit_lower_parameter_boundary=best_parameter_value
        <= lower + max((upper - lower) / 160.0, 1e-6),
        hit_upper_parameter_boundary=best_parameter_value
        >= upper - max((upper - lower) / 160.0, 1e-6),
        transformed_tree_is_ultrametric=transformed_ultrametric.ultrametric,
        transformed_tree_minimum_tip_depth=transformed_ultrametric.minimum_tip_depth,
        transformed_tree_maximum_tip_depth=transformed_ultrametric.maximum_tip_depth,
        profile_rows=[
            DiscreteMkTransformProfileRow(
                transform_parameter_value=transform_parameter_value,
                log_likelihood=log_likelihood,
            )
            for transform_parameter_value, log_likelihood in sorted(
                {row[0]: row[1] for row in profile_rows}.items()
            )
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
    transformed_tree = transform_tree_for_evolutionary_mode(
        dataset.tree,
        mode=_discrete_transform_mode_name(transform),
        parameter_value=parameter_value,
        sigsq=1.0,
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


def _discrete_transform_mode_name(transform: str) -> str:
    if transform == "lambda":
        return "pagel-lambda"
    if transform == "kappa":
        return "pagel-kappa"
    raise ComparativeMethodError(
        "unsupported discrete Mk transform; expected one of: lambda, kappa"
    )


def _transform_cache_key(parameter_value: float) -> float:
    return float(format(parameter_value, ".15g"))


def _grid_values(lower: float, upper: float, *, point_count: int) -> list[float]:
    if point_count < 2:
        return [lower]
    step = (upper - lower) / (point_count - 1)
    return [lower + (step * index) for index in range(point_count)]


def _transform_search_bracket(
    *,
    coarse_candidates: list[float],
    best_index: int,
    lower: float,
    upper: float,
) -> tuple[float, float]:
    if best_index <= 0:
        return lower, coarse_candidates[1]
    if best_index >= len(coarse_candidates) - 1:
        return coarse_candidates[-2], upper
    return coarse_candidates[best_index - 1], coarse_candidates[best_index + 1]


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
                kind=f"boundary_{transform}",
                message=f"best-supported discrete Mk {transform} falls on the search boundary and may not be well identified",
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
        pruned_missing_value_taxa=list(dataset.dropped_missing_taxa),
        warnings=warnings,
    )


def _fit_aicc(aic: float, *, sample_size: int, parameter_count: int) -> float:
    denominator = sample_size - parameter_count - 1
    if denominator <= 0:
        return math.inf
    return aic + ((2.0 * parameter_count * (parameter_count + 1.0)) / denominator)
