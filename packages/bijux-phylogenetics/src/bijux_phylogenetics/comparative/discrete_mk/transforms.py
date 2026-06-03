from __future__ import annotations

import math

from bijux_phylogenetics.ancestral.common import AncestralDiscreteDataset
from bijux_phylogenetics.ancestral.discrete import DiscreteOptimizerDiagnostics
from bijux_phylogenetics.ancestral.discrete.likelihood.fit_workflow import (
    fit_discrete_mk_model as _fit_discrete_mk_model,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    tree_log_likelihood as _tree_log_likelihood,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    parameter_count as _parameter_count,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_allowed_transition_pairs as _resolve_allowed_transition_pairs,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)
from bijux_phylogenetics.comparative.common import tip_root_depths
from bijux_phylogenetics.comparative.evolutionary_modes import (
    transform_tree_for_evolutionary_mode,
)
from bijux_phylogenetics.comparative.search import (
    BoundedSearchControls,
    run_bounded_maximization,
)
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .models import (
    DISCRETE_MK_MODEL_COMPARISON_ORDER,
    DiscreteMkTransformFit,
    DiscreteMkTransformProfileRow,
    DiscreteMkTransformWarning,
)

DISCRETE_TRANSFORM_COARSE_GRID_POINT_COUNT = 9
DISCRETE_TRANSFORM_FINE_GRID_POINT_COUNT = 17
DISCRETE_TRANSFORM_REFINEMENT_START_COUNT = 1
DISCRETE_DELTA_LOWER_BOUND = math.exp(-5.0)
DISCRETE_DELTA_UPPER_BOUND = 3.0
DISCRETE_EARLY_BURST_LOWER_BOUND = -10.0
DISCRETE_EARLY_BURST_UPPER_BOUND = 10.0


def resolve_discrete_transform_name(transform: str | None) -> str | None:
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


def validate_discrete_transform_request(
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
        validate_lambda_bounds(lambda_bounds)
        return
    if transform == "kappa":
        validate_kappa_bounds(kappa_bounds)
        return
    if transform == "delta":
        validate_delta_bounds(delta_bounds)
        return
    validate_early_burst_bounds(early_burst_bounds)


def validate_lambda_bounds(bounds: tuple[float, float]) -> None:
    lower, upper = bounds
    if not 0.0 <= lower < upper <= 1.0:
        raise ComparativeMethodError(
            "discrete Mk lambda bounds must be strictly increasing within [0, 1]"
        )


def validate_kappa_bounds(bounds: tuple[float, float]) -> None:
    lower, upper = bounds
    if not 0.0 <= lower < upper <= 1.0:
        raise ComparativeMethodError(
            "discrete Mk kappa bounds must be strictly increasing within [0, 1]"
        )


def validate_delta_bounds(bounds: tuple[float, float]) -> None:
    lower, upper = bounds
    if not 0.0 < lower < upper <= DISCRETE_DELTA_UPPER_BOUND:
        raise ComparativeMethodError(
            "discrete Mk delta bounds must be strictly increasing within "
            f"(0, {DISCRETE_DELTA_UPPER_BOUND:g}]"
        )


def validate_early_burst_bounds(bounds: tuple[float, float]) -> None:
    lower, upper = bounds
    if not (
        DISCRETE_EARLY_BURST_LOWER_BOUND
        <= lower
        < upper
        <= DISCRETE_EARLY_BURST_UPPER_BOUND
    ):
        raise ComparativeMethodError(
            "discrete Mk early-burst bounds must be strictly increasing within "
            f"[{DISCRETE_EARLY_BURST_LOWER_BOUND:g}, {DISCRETE_EARLY_BURST_UPPER_BOUND:g}]"
        )


def comparison_models(models: tuple[str, ...] | None) -> tuple[str, ...]:
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


def discrete_parameter_count(
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


def comparison_parameter_count(
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
    return discrete_parameter_count(
        state_count=len(state_order),
        model=model,
        transform=resolve_discrete_transform_name(transform),
        state_ordering=state_ordering,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
    )


def fit_discrete_mk_surface(
    dataset: AncestralDiscreteDataset,
    *,
    model: str,
    ascertainment_policy: str,
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
            ascertainment_policy=ascertainment_policy,
        )
        return dataset.tree, rate_matrix, root_prior, optimizer_diagnostics, None, []
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
    return fit_discrete_mk_parameterized_transform_surface(
        dataset,
        model=model,
        ascertainment_policy=ascertainment_policy,
        transform=transform,
        state_ordering=state_ordering,
        state_order=state_order,
        allowed_transition_pairs=allowed_transition_pairs,
        bounds=bounds,
    )


def fit_discrete_mk_parameterized_transform_surface(
    dataset: AncestralDiscreteDataset,
    *,
    model: str,
    ascertainment_policy: str,
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
            coarse_grid_point_count=DISCRETE_TRANSFORM_COARSE_GRID_POINT_COUNT,
            fine_grid_point_count=DISCRETE_TRANSFORM_FINE_GRID_POINT_COUNT,
            refinement_start_count=DISCRETE_TRANSFORM_REFINEMENT_START_COUNT,
        ),
        evaluate=lambda parameter_value: (
            evaluate_discrete_mk_transform_search_candidate(
                dataset,
                model=model,
                ascertainment_policy=ascertainment_policy,
                transform=transform,
                state_ordering=state_ordering,
                state_order=state_order,
                allowed_transition_pairs=allowed_transition_pairs,
                parameter_value=parameter_value,
            )
        ),
        optimizer_name="governed-two-stage-grid-search",
        parameter_search_strategy="bounded-two-stage-grid-search",
    )
    transformed_tree, rate_matrix, root_prior, optimizer_diagnostics = (
        search_result.payload
    )
    best_parameter_value = search_result.parameter_value
    transform_warning_rows = discrete_transform_warning_rows(
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
        parameter_name=discrete_transform_parameter_name(transform),
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


def evaluate_discrete_mk_transform_candidate(
    dataset: AncestralDiscreteDataset,
    *,
    model: str,
    ascertainment_policy: str,
    transform: str,
    state_ordering: str,
    state_order: list[str],
    allowed_transition_pairs: set[tuple[int, int]],
    parameter_value: float,
) -> tuple[object, object, object, DiscreteOptimizerDiagnostics, float]:
    parameter_value = transform_cache_key(parameter_value)
    transformed_tree = transform_discrete_mk_tree(
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
        ascertainment_policy=ascertainment_policy,
    )
    log_likelihood = _tree_log_likelihood(
        transformed_tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        root_prior_mode="observed",
        ascertainment_policy=ascertainment_policy,
    )
    return (
        transformed_tree,
        rate_matrix,
        root_prior,
        optimizer_diagnostics,
        log_likelihood,
    )


def evaluate_discrete_mk_transform_search_candidate(
    dataset: AncestralDiscreteDataset,
    *,
    model: str,
    ascertainment_policy: str,
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
    ) = evaluate_discrete_mk_transform_candidate(
        dataset,
        model=model,
        ascertainment_policy=ascertainment_policy,
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


def discrete_transform_mode_name(transform: str) -> str:
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


def discrete_transform_parameter_name(transform: str) -> str:
    if transform == "early-burst":
        return "a"
    return transform


def transform_discrete_mk_tree(
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
        mode=discrete_transform_mode_name(transform),
        parameter_value=transformed_parameter_value,
        sigsq=1.0,
    )


def transform_cache_key(parameter_value: float) -> float:
    return float(format(parameter_value, ".15g"))


def discrete_transform_warning_rows(
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
    if math.isclose(
        parameter_value, lower, rel_tol=0.0, abs_tol=boundary_tolerance
    ) or math.isclose(parameter_value, upper, rel_tol=0.0, abs_tol=boundary_tolerance):
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
    if transform == "lambda" and parameter_value <= lower + max(
        boundary_tolerance, 1e-6
    ):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="weak_phylogenetic_signal",
                message="best-supported discrete Mk lambda remains close to the zero-signal boundary and may be difficult to distinguish from a star-like transition surface",
            )
        )
    if transform == "lambda" and parameter_value >= upper - max(
        boundary_tolerance, 1e-6
    ):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="brownian_limit",
                message="best-supported discrete Mk lambda remains close to the untransformed boundary and may be difficult to distinguish from the original branch-length surface",
            )
        )
    if transform == "kappa" and parameter_value <= lower + max(
        boundary_tolerance, 1e-6
    ):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="branch_length_flattening_limit",
                message="best-supported discrete Mk kappa remains close to the zero-contrast boundary and may be difficult to distinguish from branch-length flattening",
            )
        )
    if transform == "kappa" and parameter_value >= upper - max(
        boundary_tolerance, 1e-6
    ):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="branch_length_contrast_limit",
                message="best-supported discrete Mk kappa remains close to the untransformed branch-length boundary and may be weakly identified",
            )
        )
    if transform == "delta" and parameter_value <= lower + max(
        boundary_tolerance, 1e-6
    ):
        warnings.append(
            DiscreteMkTransformWarning(
                kind="earliest_change_limit",
                message="best-supported discrete Mk delta remains close to the earliest-change boundary and may be difficult to distinguish from an extreme root-concentrated transition surface",
            )
        )
    if transform == "delta" and parameter_value >= upper - max(
        boundary_tolerance, 1e-6
    ):
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
