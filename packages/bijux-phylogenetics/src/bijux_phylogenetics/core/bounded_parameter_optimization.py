from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Generic, TypeVar

from .bounded_parameter_validation import (
    validate_increasing_parameter_bounds,
    validate_parameter_within_bounds,
)
from .bounded_scalar_search import (
    BoundedScalarSearchControls,
    BoundedScalarSearchDiagnostics,
    BoundedScalarSearchProfileRow,
    run_bounded_scalar_golden_section_maximization,
    run_bounded_scalar_maximization,
)

_PayloadT = TypeVar("_PayloadT")

BoundedParameterSearchControls = BoundedScalarSearchControls
BoundedParameterSearchDiagnostics = BoundedScalarSearchDiagnostics
BoundedParameterSearchProfileRow = BoundedScalarSearchProfileRow


@dataclass(slots=True)
class BoundedParameterOptimizationResult(Generic[_PayloadT]):
    """Best result from one shared bounded one-parameter optimization."""

    parameter_value: float
    payload: _PayloadT
    objective_value: float
    function_evaluation_count: int
    converged: bool
    diagnostics: BoundedParameterSearchDiagnostics
    profile_rows: list[BoundedParameterSearchProfileRow]


@dataclass(slots=True)
class BoundedCoordinateOptimizationResult(Generic[_PayloadT]):
    """Best result from one shared bounded coordinate-wise optimization."""

    parameter_values: dict[str, float]
    payload: _PayloadT
    objective_value: float
    function_evaluation_count: int
    optimization_pass_count: int
    converged: bool
    parameter_diagnostics: dict[str, BoundedParameterSearchDiagnostics]


def run_bounded_parameter_maximization(
    *,
    lower_bound: float,
    upper_bound: float,
    controls: BoundedParameterSearchControls,
    evaluate: Callable[[float], tuple[_PayloadT, float]],
    parameter_name: str = "search parameter",
    owner_name: str = "bounded parameter search",
    optimizer_name: str = "governed-bounded-multi-start-search",
    parameter_search_strategy: str = "bounded-coarse-grid-with-multi-start-golden-section-refinement",
) -> BoundedParameterOptimizationResult[_PayloadT]:
    """Maximize one bounded scalar objective through the shared optimizer owner."""
    validated_lower_bound, validated_upper_bound = validate_increasing_parameter_bounds(
        parameter_name=parameter_name,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        owner_name=owner_name,
    )

    def evaluate_validated(parameter_value: float) -> tuple[_PayloadT, float]:
        validated_parameter_value = validate_parameter_within_bounds(
            parameter_name=parameter_name,
            value=parameter_value,
            lower_bound=validated_lower_bound,
            upper_bound=validated_upper_bound,
            owner_name=owner_name,
        )
        return evaluate(validated_parameter_value)

    search_result = run_bounded_scalar_maximization(
        lower_bound=validated_lower_bound,
        upper_bound=validated_upper_bound,
        controls=controls,
        evaluate=evaluate_validated,
        optimizer_name=optimizer_name,
        parameter_search_strategy=parameter_search_strategy,
    )
    return BoundedParameterOptimizationResult(
        parameter_value=search_result.parameter_value,
        payload=search_result.payload,
        objective_value=search_result.objective_value,
        function_evaluation_count=search_result.diagnostics.function_evaluation_count,
        converged=search_result.diagnostics.converged,
        diagnostics=search_result.diagnostics,
        profile_rows=search_result.profile_rows,
    )


def run_bounded_parameter_golden_section_maximization(
    *,
    lower_bound: float,
    upper_bound: float,
    evaluate: Callable[[float], tuple[_PayloadT, float]],
    parameter_name: str = "search parameter",
    owner_name: str = "bounded parameter search",
    tolerance: float = 1e-9,
    max_iterations: int = 400,
    optimizer_name: str = "golden-section-search",
    parameter_search_strategy: str = "bounded-single-start-golden-section-search",
) -> BoundedParameterOptimizationResult[_PayloadT]:
    """Maximize one bounded scalar objective by shared golden-section search."""
    validated_lower_bound, validated_upper_bound = validate_increasing_parameter_bounds(
        parameter_name=parameter_name,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        owner_name=owner_name,
    )

    def evaluate_validated(parameter_value: float) -> tuple[_PayloadT, float]:
        validated_parameter_value = validate_parameter_within_bounds(
            parameter_name=parameter_name,
            value=parameter_value,
            lower_bound=validated_lower_bound,
            upper_bound=validated_upper_bound,
            owner_name=owner_name,
        )
        return evaluate(validated_parameter_value)

    search_result = run_bounded_scalar_golden_section_maximization(
        lower_bound=validated_lower_bound,
        upper_bound=validated_upper_bound,
        evaluate=evaluate_validated,
        tolerance=tolerance,
        max_iterations=max_iterations,
        optimizer_name=optimizer_name,
        parameter_search_strategy=parameter_search_strategy,
    )
    return BoundedParameterOptimizationResult(
        parameter_value=search_result.parameter_value,
        payload=search_result.payload,
        objective_value=search_result.objective_value,
        function_evaluation_count=search_result.diagnostics.function_evaluation_count,
        converged=search_result.diagnostics.converged,
        diagnostics=search_result.diagnostics,
        profile_rows=search_result.profile_rows,
    )


def run_bounded_coordinate_optimization(
    *,
    initial_values: dict[str, float],
    bounds_by_name: dict[str, tuple[float, float]],
    evaluate: Callable[[dict[str, float]], tuple[_PayloadT, float]],
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
    owner_name: str = "coordinate parameter search",
) -> BoundedCoordinateOptimizationResult[_PayloadT]:
    """Maximize one bounded coordinate-wise objective through the shared owner."""
    if max_coordinate_passes < 1:
        raise ValueError(f"{owner_name} requires at least one optimization pass")
    if improvement_tolerance < 0.0:
        raise ValueError(f"{owner_name} requires a nonnegative improvement tolerance")
    if set(initial_values) != set(bounds_by_name):
        raise ValueError(
            f"{owner_name} requires identical parameter names in starts and bounds"
        )

    current_values = {name: float(value) for name, value in initial_values.items()}
    validated_bounds_by_name: dict[str, tuple[float, float]] = {}
    for name, (lower_bound, upper_bound) in bounds_by_name.items():
        validated_lower_bound, validated_upper_bound = (
            validate_increasing_parameter_bounds(
                parameter_name=name,
                lower_bound=lower_bound,
                upper_bound=upper_bound,
                owner_name=owner_name,
            )
        )
        validated_bounds_by_name[name] = (
            validated_lower_bound,
            validated_upper_bound,
        )
        current_values[name] = validate_parameter_within_bounds(
            parameter_name=name,
            value=current_values[name],
            lower_bound=validated_lower_bound,
            upper_bound=validated_upper_bound,
            owner_name=owner_name,
        )

    current_payload, current_objective = evaluate(dict(current_values))
    function_evaluation_count = 1
    optimization_pass_count = 0
    converged = False
    all_coordinate_searches_converged = True
    parameter_diagnostics: dict[str, BoundedParameterSearchDiagnostics] = {}

    for optimization_pass in range(1, max_coordinate_passes + 1):
        optimization_pass_count = optimization_pass
        improved = False
        for name in current_values:
            lower_bound, upper_bound = validated_bounds_by_name[name]

            def evaluate_coordinate(
                candidate_value: float,
                name: str = name,
            ) -> tuple[_PayloadT, float]:
                candidate_values = dict(current_values)
                candidate_values[name] = candidate_value
                return evaluate(candidate_values)

            search_result = run_bounded_parameter_golden_section_maximization(
                lower_bound=lower_bound,
                upper_bound=upper_bound,
                evaluate=evaluate_coordinate,
                parameter_name=name,
                owner_name=owner_name,
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

    return BoundedCoordinateOptimizationResult(
        parameter_values=dict(current_values),
        payload=current_payload,
        objective_value=current_objective,
        function_evaluation_count=function_evaluation_count,
        optimization_pass_count=optimization_pass_count,
        converged=converged,
        parameter_diagnostics=parameter_diagnostics,
    )


__all__ = [
    "BoundedCoordinateOptimizationResult",
    "BoundedParameterOptimizationResult",
    "BoundedParameterSearchControls",
    "BoundedParameterSearchDiagnostics",
    "BoundedParameterSearchProfileRow",
    "run_bounded_coordinate_optimization",
    "run_bounded_parameter_golden_section_maximization",
    "run_bounded_parameter_maximization",
]
