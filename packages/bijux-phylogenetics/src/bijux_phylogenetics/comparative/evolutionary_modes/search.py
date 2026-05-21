from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    build_brownian_covariance_matrix,
    stable_covariance,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    _fit_intercept_only_model,
)
from bijux_phylogenetics.comparative.evolutionary_modes.models import (
    ContinuousModeBoundaryAssessment,
    ContinuousModeOptimizerDiagnostics,
    ContinuousModeSearchControls,
    EvolutionaryModeIdentifiabilityWarning,
)
from bijux_phylogenetics.comparative.evolutionary_modes.numeric import stable_float
from bijux_phylogenetics.comparative.evolutionary_modes.tree_transforms import (
    transform_tree,
)
from bijux_phylogenetics.comparative.search import (
    BoundedSearchControls,
    run_bounded_maximization,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


@dataclass(slots=True)
class TransformedModeSearchResult:
    parameter_value: float
    transformed_tree: PhyloTree
    covariance: list[list[float]]
    optimizer_diagnostics: ContinuousModeOptimizerDiagnostics
    profile: list[tuple[float, float]]


def reject_nonparameterized_search_controls(
    mode: str,
    search_controls: ContinuousModeSearchControls | None,
) -> None:
    if search_controls is None:
        return
    raise ComparativeMethodError(
        f"{mode} mode does not expose bounded parameter-search controls because it uses a closed-form profile solution"
    )


def normalized_search_controls(
    search_controls: ContinuousModeSearchControls | None,
) -> ContinuousModeSearchControls:
    controls = (
        ContinuousModeSearchControls() if search_controls is None else search_controls
    )
    if controls.coarse_grid_point_count < 2:
        raise ComparativeMethodError(
            "coarse_grid_point_count must be at least 2 for bounded parameter search"
        )
    if controls.fine_grid_point_count < 2:
        raise ComparativeMethodError(
            "fine_grid_point_count must be at least 2 for bounded parameter search"
        )
    if controls.refinement_start_count < 1:
        raise ComparativeMethodError(
            "refinement_start_count must be at least 1 for bounded parameter search"
        )
    return controls


def best_transformed_mode_fit(
    dataset: ComparativeDataset,
    *,
    mode: str,
    bounds: tuple[float, float],
    search_controls: ContinuousModeSearchControls,
) -> TransformedModeSearchResult:
    lower, upper = bounds
    if upper <= lower:
        raise ComparativeMethodError("parameter bounds must be strictly increasing")
    if mode == "ornstein-uhlenbeck":
        lower = max(lower, 1e-6)
    search_result = run_bounded_maximization(
        lower_bound=lower,
        upper_bound=upper,
        controls=BoundedSearchControls(
            coarse_grid_point_count=search_controls.coarse_grid_point_count,
            fine_grid_point_count=search_controls.fine_grid_point_count,
            initial_parameter_value=search_controls.initial_parameter_value,
            refinement_start_count=search_controls.refinement_start_count,
        ),
        evaluate=lambda parameter_value: evaluate_transformed_mode_candidate(
            dataset,
            mode=mode,
            parameter_value=parameter_value,
        ),
        optimizer_name="governed-multi-start-grid-search",
        parameter_search_strategy="bounded-coarse-grid-with-multi-start-golden-section-refinement",
    )
    best_tree, best_covariance = search_result.payload
    diagnostics = ContinuousModeOptimizerDiagnostics(
        optimizer_name="governed-two-stage-grid-search",
        parameter_search_strategy="bounded-two-stage-grid-search",
        lower_bound=stable_float(search_result.diagnostics.lower_bound),
        upper_bound=stable_float(search_result.diagnostics.upper_bound),
        starting_parameter_policy=search_result.diagnostics.starting_parameter_policy,
        starting_parameter_value=stable_float(
            search_result.diagnostics.starting_parameter_value
        ),
        starting_parameter_log_likelihood=stable_float(
            search_result.diagnostics.starting_parameter_objective_value
        ),
        coarse_grid_point_count=search_result.diagnostics.coarse_grid_point_count,
        fine_grid_point_count=search_result.diagnostics.fine_grid_point_count,
        refinement_start_count=search_result.diagnostics.refinement_start_count,
        function_evaluation_count=search_result.diagnostics.function_evaluation_count,
        coarse_best_parameter=stable_float(
            search_result.diagnostics.coarse_best_parameter
        ),
        coarse_best_log_likelihood=stable_float(
            search_result.diagnostics.coarse_best_objective_value
        ),
        fine_search_start=stable_float(search_result.diagnostics.fine_search_start),
        fine_search_stop=stable_float(search_result.diagnostics.fine_search_stop),
        converged=search_result.diagnostics.converged,
        hit_lower_boundary=search_result.diagnostics.hit_lower_boundary,
        hit_upper_boundary=search_result.diagnostics.hit_upper_boundary,
    )
    return TransformedModeSearchResult(
        parameter_value=search_result.parameter_value,
        transformed_tree=best_tree,
        covariance=best_covariance,
        optimizer_diagnostics=diagnostics,
        profile=[
            (row.parameter_value, row.objective_value)
            for row in search_result.profile_rows
        ],
    )


def best_pagel_lambda_fit(
    dataset: ComparativeDataset,
    *,
    bounds: tuple[float, float],
    search_controls: ContinuousModeSearchControls,
) -> TransformedModeSearchResult:
    lower, upper = bounds
    if lower < 0.0 or upper > 1.0 or upper <= lower:
        raise ComparativeMethodError(
            "Pagel-lambda bounds must be strictly increasing within [0, 1]"
        )
    search_result = run_bounded_maximization(
        lower_bound=lower,
        upper_bound=upper,
        controls=BoundedSearchControls(
            coarse_grid_point_count=search_controls.coarse_grid_point_count,
            fine_grid_point_count=search_controls.fine_grid_point_count,
            initial_parameter_value=search_controls.initial_parameter_value,
            refinement_start_count=search_controls.refinement_start_count,
        ),
        evaluate=lambda parameter_value: evaluate_transformed_mode_candidate(
            dataset,
            mode="pagel-lambda",
            parameter_value=parameter_value,
        ),
        optimizer_name="governed-multi-start-grid-search",
        parameter_search_strategy="bounded-coarse-grid-with-multi-start-golden-section-refinement",
    )
    best_tree, best_covariance = search_result.payload
    diagnostics = ContinuousModeOptimizerDiagnostics(
        optimizer_name="governed-two-stage-grid-search",
        parameter_search_strategy="bounded-two-stage-grid-search",
        lower_bound=stable_float(search_result.diagnostics.lower_bound),
        upper_bound=stable_float(search_result.diagnostics.upper_bound),
        starting_parameter_policy=search_result.diagnostics.starting_parameter_policy,
        starting_parameter_value=stable_float(
            search_result.diagnostics.starting_parameter_value
        ),
        starting_parameter_log_likelihood=stable_float(
            search_result.diagnostics.starting_parameter_objective_value
        ),
        coarse_grid_point_count=search_result.diagnostics.coarse_grid_point_count,
        fine_grid_point_count=search_result.diagnostics.fine_grid_point_count,
        refinement_start_count=search_result.diagnostics.refinement_start_count,
        function_evaluation_count=search_result.diagnostics.function_evaluation_count,
        coarse_best_parameter=stable_float(
            search_result.diagnostics.coarse_best_parameter
        ),
        coarse_best_log_likelihood=stable_float(
            search_result.diagnostics.coarse_best_objective_value
        ),
        fine_search_start=stable_float(search_result.diagnostics.fine_search_start),
        fine_search_stop=stable_float(search_result.diagnostics.fine_search_stop),
        converged=search_result.diagnostics.converged,
        hit_lower_boundary=search_result.diagnostics.hit_lower_boundary,
        hit_upper_boundary=search_result.diagnostics.hit_upper_boundary,
    )
    return TransformedModeSearchResult(
        parameter_value=search_result.parameter_value,
        transformed_tree=best_tree,
        covariance=best_covariance,
        optimizer_diagnostics=diagnostics,
        profile=[
            (row.parameter_value, row.objective_value)
            for row in search_result.profile_rows
        ],
    )


def evaluate_transformed_mode_candidate(
    dataset: ComparativeDataset,
    *,
    mode: str,
    parameter_value: float,
) -> tuple[tuple[PhyloTree, list[list[float]]], float]:
    transformed_tree = transform_tree(
        dataset.tree,
        mode=mode,
        parameter_value=parameter_value,
    )
    covariance = stable_covariance(
        build_brownian_covariance_matrix(transformed_tree, dataset.taxa)
    )
    fit = _fit_intercept_only_model(dataset, covariance)
    return (transformed_tree, covariance), fit.log_likelihood


def ou_identifiability_warnings_from_profile(
    dataset: ComparativeDataset,
    alpha: float,
    profile: list[tuple[float, float]],
) -> list[EvolutionaryModeIdentifiabilityWarning]:
    warnings: list[EvolutionaryModeIdentifiabilityWarning] = []
    if len(dataset.taxa) < 5:
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="small_sample_size",
                message="OU alpha is hard to identify with fewer than five taxa",
            )
        )
    ordered_alphas = sorted(candidate for candidate, _ in profile)
    if math.isclose(
        alpha, ordered_alphas[0], rel_tol=0.0, abs_tol=1e-9
    ) or math.isclose(alpha, ordered_alphas[-1], rel_tol=0.0, abs_tol=1e-9):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="boundary_alpha",
                message="best-supported OU alpha falls on the search boundary and may not be well identified",
            )
        )
    ordered_log_likelihoods = sorted(
        (log_likelihood for _, log_likelihood in profile),
        reverse=True,
    )
    if len(ordered_log_likelihoods) > 1 and (
        ordered_log_likelihoods[0] - ordered_log_likelihoods[1] < 0.5
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="flat_likelihood",
                message="OU likelihood surface is shallow across alpha values, so model choice may be unstable",
            )
        )
    alpha_span = ordered_alphas[-1] - ordered_alphas[0]
    weak_pull_threshold = ordered_alphas[0] + (alpha_span / 3.0)
    if alpha <= weak_pull_threshold:
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="weak_pull_to_optimum",
                message="best-supported OU alpha is weak and may be difficult to distinguish from Brownian motion",
            )
        )
    return warnings


def early_burst_identifiability_warnings_from_profile(
    rate_change: float,
    profile: list[tuple[float, float]],
    bounds: tuple[float, float],
) -> list[EvolutionaryModeIdentifiabilityWarning]:
    lower, upper = bounds
    span = upper - lower
    ordered_log_likelihoods = sorted(
        (log_likelihood for _, log_likelihood in profile),
        reverse=True,
    )
    warnings: list[EvolutionaryModeIdentifiabilityWarning] = []
    boundary_tolerance = max(span / 160.0, 1e-9)
    if math.isclose(
        rate_change, lower, rel_tol=0.0, abs_tol=boundary_tolerance
    ) or math.isclose(rate_change, upper, rel_tol=0.0, abs_tol=boundary_tolerance):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="boundary_rate_change",
                message="best-supported early-burst rate change falls on the search boundary and may not be well identified",
            )
        )
    if len(ordered_log_likelihoods) > 1 and (
        ordered_log_likelihoods[0] - ordered_log_likelihoods[1] < 0.5
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="flat_likelihood_profile",
                message="early-burst likelihood stays shallow across rate-change values, so model choice may be unstable",
            )
        )
    if rate_change <= lower + max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="brownian_like_rate_change",
                message="best-supported early-burst rate change remains close to the zero-change boundary and may be difficult to distinguish from Brownian motion",
            )
        )
    return warnings


def lambda_identifiability_warnings_from_profile(
    lambda_value: float,
    profile: list[tuple[float, float]],
    bounds: tuple[float, float],
) -> list[EvolutionaryModeIdentifiabilityWarning]:
    lower, upper = bounds
    span = upper - lower
    ordered_log_likelihoods = sorted(
        (log_likelihood for _, log_likelihood in profile),
        reverse=True,
    )
    warnings: list[EvolutionaryModeIdentifiabilityWarning] = []
    boundary_tolerance = max(span / 160.0, 1e-9)
    if math.isclose(
        lambda_value, lower, rel_tol=0.0, abs_tol=boundary_tolerance
    ) or math.isclose(lambda_value, upper, rel_tol=0.0, abs_tol=boundary_tolerance):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="boundary_lambda",
                message="best-supported Pagel lambda falls on the search boundary and may not be well identified",
            )
        )
    if len(ordered_log_likelihoods) > 1 and (
        ordered_log_likelihoods[0] - ordered_log_likelihoods[1] < 0.5
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="flat_likelihood",
                message="Pagel-lambda likelihood stays shallow across the bounded search, so the covariance scaling may be unstable",
            )
        )
    if lambda_value <= lower + max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="weak_phylogenetic_signal",
                message="best-supported Pagel lambda remains close to the zero-signal boundary and may be difficult to distinguish from a star-like covariance surface",
            )
        )
    if lambda_value >= upper - max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="brownian_limit",
                message="best-supported Pagel lambda remains close to the Brownian boundary and may be difficult to distinguish from an untransformed covariance surface",
            )
        )
    return warnings


def kappa_identifiability_warnings_from_profile(
    kappa_value: float,
    profile: list[tuple[float, float]],
    bounds: tuple[float, float],
) -> list[EvolutionaryModeIdentifiabilityWarning]:
    lower, upper = bounds
    span = upper - lower
    ordered_log_likelihoods = sorted(
        (log_likelihood for _, log_likelihood in profile),
        reverse=True,
    )
    warnings: list[EvolutionaryModeIdentifiabilityWarning] = []
    boundary_tolerance = max(span / 160.0, 1e-9)
    if math.isclose(
        kappa_value, lower, rel_tol=0.0, abs_tol=boundary_tolerance
    ) or math.isclose(kappa_value, upper, rel_tol=0.0, abs_tol=boundary_tolerance):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="boundary_kappa",
                message="best-supported Pagel kappa falls on the search boundary and may not be well identified",
            )
        )
    if len(ordered_log_likelihoods) > 1 and (
        ordered_log_likelihoods[0] - ordered_log_likelihoods[1] < 0.5
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="flat_likelihood",
                message="Pagel-kappa likelihood stays shallow across the bounded search, so branch-length transformation support may be unstable",
            )
        )
    if kappa_value <= lower + max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="punctuational_limit",
                message="best-supported Pagel kappa remains close to the equal-length punctuational boundary and may be difficult to distinguish from a branch-count surface",
            )
        )
    if kappa_value >= upper - max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="upper_search_limit",
                message="best-supported Pagel kappa remains close to the upper search boundary, so a wider branch-length review may be needed",
            )
        )
    return warnings


def delta_identifiability_warnings_from_profile(
    delta_value: float,
    profile: list[tuple[float, float]],
    bounds: tuple[float, float],
) -> list[EvolutionaryModeIdentifiabilityWarning]:
    lower, upper = bounds
    span = upper - lower
    ordered_log_likelihoods = sorted(
        (log_likelihood for _, log_likelihood in profile),
        reverse=True,
    )
    warnings: list[EvolutionaryModeIdentifiabilityWarning] = []
    boundary_tolerance = max(span / 160.0, 1e-9)
    if math.isclose(
        delta_value, lower, rel_tol=0.0, abs_tol=boundary_tolerance
    ) or math.isclose(delta_value, upper, rel_tol=0.0, abs_tol=boundary_tolerance):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="boundary_delta",
                message="best-supported Pagel delta falls on the search boundary and may not be well identified",
            )
        )
    if len(ordered_log_likelihoods) > 1 and (
        ordered_log_likelihoods[0] - ordered_log_likelihoods[1] < 0.5
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="flat_likelihood",
                message="Pagel-delta likelihood stays shallow across the bounded search, so temporal concentration support may be unstable",
            )
        )
    if delta_value <= lower + max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="early_change_limit",
                message="best-supported Pagel delta remains close to the earliest-change boundary and may be difficult to distinguish from an extreme root-concentrated surface",
            )
        )
    if delta_value >= upper - max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="late_change_limit",
                message="best-supported Pagel delta remains close to the late-change boundary and may be difficult to distinguish from an extreme tip-concentrated surface",
            )
        )
    return warnings


def continuous_boundary_assessment(
    *,
    parameter_name: str | None,
    parameter_value: float | None,
    optimizer_diagnostics: ContinuousModeOptimizerDiagnostics | None,
    identifiability_warnings: list[EvolutionaryModeIdentifiabilityWarning],
) -> ContinuousModeBoundaryAssessment | None:
    if (
        parameter_name is None
        or parameter_value is None
        or optimizer_diagnostics is None
    ):
        return None
    lower_bound = optimizer_diagnostics.lower_bound
    upper_bound = optimizer_diagnostics.upper_bound
    near_boundary_tolerance = max((upper_bound - lower_bound) / 20.0, 1e-6)
    warning_kinds = [warning.kind for warning in identifiability_warnings]
    boundary_warning_kinds = [
        kind for kind in warning_kinds if is_boundary_warning_kind(kind)
    ]
    hit_lower_boundary = optimizer_diagnostics.hit_lower_boundary
    hit_upper_boundary = optimizer_diagnostics.hit_upper_boundary
    near_lower_boundary = hit_lower_boundary or (
        parameter_value <= lower_bound + near_boundary_tolerance
    )
    near_upper_boundary = hit_upper_boundary or (
        parameter_value >= upper_bound - near_boundary_tolerance
    )
    flat_likelihood_near_boundary = any(
        kind in {"flat_likelihood", "flat_likelihood_profile"}
        for kind in boundary_warning_kinds
    )
    boundary_dominates_interpretation = (
        near_lower_boundary or near_upper_boundary
    ) and (
        flat_likelihood_near_boundary
        or any(
            kind.startswith("boundary_") or kind in boundary_limit_warning_kinds()
            for kind in boundary_warning_kinds
        )
    )
    return ContinuousModeBoundaryAssessment(
        affected_parameter=parameter_name,
        parameter_value=stable_float(parameter_value),
        lower_bound=stable_float(lower_bound),
        upper_bound=stable_float(upper_bound),
        hit_lower_boundary=hit_lower_boundary,
        hit_upper_boundary=hit_upper_boundary,
        near_lower_boundary=near_lower_boundary,
        near_upper_boundary=near_upper_boundary,
        flat_likelihood_near_boundary=flat_likelihood_near_boundary,
        boundary_warning_kinds=boundary_warning_kinds,
        boundary_dominates_interpretation=boundary_dominates_interpretation,
        stable_conclusion_supported=not boundary_dominates_interpretation,
    )


def is_boundary_warning_kind(kind: str) -> bool:
    return kind.startswith("boundary_") or kind in {
        "flat_likelihood",
        "flat_likelihood_profile",
        *boundary_limit_warning_kinds(),
    }


def boundary_limit_warning_kinds() -> set[str]:
    return {
        "weak_pull_to_optimum",
        "brownian_like_rate_change",
        "weak_phylogenetic_signal",
        "brownian_limit",
        "punctuational_limit",
        "upper_search_limit",
        "early_change_limit",
        "late_change_limit",
    }
