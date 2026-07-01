from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Generic, TypeVar

from bijux_phylogenetics.core.bounded_scalar_search import (
    BoundedScalarSearchDiagnostics,
    run_bounded_scalar_golden_section_maximization,
)
from bijux_phylogenetics.phylo.likelihood.parameter_bounds import (
    validate_increasing_parameter_bounds,
    validate_parameter_within_bounds,
)

_PayloadT = TypeVar("_PayloadT")


@dataclass(slots=True)
class BoundedLikelihoodSearchResult(Generic[_PayloadT]):
    """Best result from one bounded one-parameter likelihood search."""

    parameter_value: float
    payload: _PayloadT
    objective_value: float
    function_evaluation_count: int
    converged: bool
    diagnostics: BoundedScalarSearchDiagnostics


@dataclass(slots=True)
class BoundedCoordinateLikelihoodSearchResult(Generic[_PayloadT]):
    """Best result from one bounded coordinate-wise likelihood search."""

    parameter_values: dict[str, float]
    payload: _PayloadT
    objective_value: float
    function_evaluation_count: int
    optimization_pass_count: int
    converged: bool
    parameter_diagnostics: dict[str, BoundedScalarSearchDiagnostics]


def run_bounded_likelihood_search(
    *,
    lower_bound: float,
    upper_bound: float,
    evaluate: Callable[[float], tuple[_PayloadT, float]],
    tolerance: float = 1e-9,
    max_iterations: int = 400,
) -> BoundedLikelihoodSearchResult[_PayloadT]:
    """Maximize one bounded likelihood-like objective by shared scalar search."""
    validated_lower_bound, validated_upper_bound = validate_increasing_parameter_bounds(
        parameter_name="search parameter",
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        owner_name="bounded likelihood search",
    )

    def evaluate_validated(parameter_value: float) -> tuple[_PayloadT, float]:
        validated_parameter_value = validate_parameter_within_bounds(
            parameter_name="search parameter",
            value=parameter_value,
            lower_bound=validated_lower_bound,
            upper_bound=validated_upper_bound,
            owner_name="bounded likelihood search",
        )
        return evaluate(validated_parameter_value)

    search_result = run_bounded_scalar_golden_section_maximization(
        lower_bound=validated_lower_bound,
        upper_bound=validated_upper_bound,
        evaluate=evaluate_validated,
        tolerance=tolerance,
        max_iterations=max_iterations,
    )
    return BoundedLikelihoodSearchResult(
        parameter_value=search_result.parameter_value,
        payload=search_result.payload,
        objective_value=search_result.objective_value,
        function_evaluation_count=search_result.diagnostics.function_evaluation_count,
        converged=search_result.diagnostics.converged,
        diagnostics=search_result.diagnostics,
    )


def run_bounded_coordinate_likelihood_search(
    *,
    initial_values: dict[str, float],
    bounds_by_name: dict[str, tuple[float, float]],
    evaluate: Callable[[dict[str, float]], tuple[_PayloadT, float]],
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> BoundedCoordinateLikelihoodSearchResult[_PayloadT]:
    """Maximize one bounded likelihood-like objective by coordinate search."""
    if max_coordinate_passes < 1:
        raise ValueError("coordinate likelihood search requires at least one pass")
    if improvement_tolerance < 0.0:
        raise ValueError(
            "coordinate likelihood search improvement_tolerance must be nonnegative"
        )
    if set(initial_values) != set(bounds_by_name):
        raise ValueError(
            "initial_values and bounds_by_name must have identical parameter names"
        )

    current_values = {name: float(value) for name, value in initial_values.items()}
    for name, (lower_bound, upper_bound) in bounds_by_name.items():
        validated_lower_bound, validated_upper_bound = (
            validate_increasing_parameter_bounds(
                parameter_name=name,
                lower_bound=lower_bound,
                upper_bound=upper_bound,
                owner_name="coordinate likelihood search",
            )
        )
        current_values[name] = validate_parameter_within_bounds(
            parameter_name=name,
            value=current_values[name],
            lower_bound=validated_lower_bound,
            upper_bound=validated_upper_bound,
            owner_name="coordinate likelihood search",
        )

    current_payload, current_objective = evaluate(dict(current_values))
    function_evaluation_count = 1
    optimization_pass_count = 0
    converged = False
    all_coordinate_searches_converged = True
    parameter_diagnostics: dict[str, BoundedScalarSearchDiagnostics] = {}

    for optimization_pass in range(1, max_coordinate_passes + 1):
        optimization_pass_count = optimization_pass
        improved = False
        for name in current_values:
            lower_bound, upper_bound = bounds_by_name[name]

            def evaluate_coordinate(
                candidate_value: float,
                name: str = name,
            ) -> tuple[_PayloadT, float]:
                candidate_values = dict(current_values)
                candidate_values[name] = candidate_value
                return evaluate(candidate_values)

            search_result = run_bounded_likelihood_search(
                lower_bound=lower_bound,
                upper_bound=upper_bound,
                evaluate=evaluate_coordinate,
            )
            parameter_diagnostics[name] = search_result.diagnostics
            function_evaluation_count += search_result.function_evaluation_count
            if not search_result.converged:
                all_coordinate_searches_converged = False
            if search_result.objective_value > (
                current_objective + improvement_tolerance
            ):
                current_values[name] = search_result.parameter_value
                current_payload = search_result.payload
                current_objective = search_result.objective_value
                improved = True
        if not improved:
            converged = all_coordinate_searches_converged
            break

    return BoundedCoordinateLikelihoodSearchResult(
        parameter_values=dict(current_values),
        payload=current_payload,
        objective_value=current_objective,
        function_evaluation_count=function_evaluation_count,
        optimization_pass_count=optimization_pass_count,
        converged=converged,
        parameter_diagnostics=parameter_diagnostics,
    )
