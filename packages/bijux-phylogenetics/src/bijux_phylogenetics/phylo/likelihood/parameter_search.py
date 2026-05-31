from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math
from typing import Generic, TypeVar

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


@dataclass(slots=True)
class BoundedCoordinateLikelihoodSearchResult(Generic[_PayloadT]):
    """Best result from one bounded coordinate-wise likelihood search."""

    parameter_values: dict[str, float]
    payload: _PayloadT
    objective_value: float
    function_evaluation_count: int
    optimization_pass_count: int
    converged: bool


def run_bounded_likelihood_search(
    *,
    lower_bound: float,
    upper_bound: float,
    evaluate: Callable[[float], tuple[_PayloadT, float]],
    tolerance: float = 1e-9,
    max_iterations: int = 400,
) -> BoundedLikelihoodSearchResult[_PayloadT]:
    """Maximize one bounded likelihood-like objective by golden-section search."""
    validated_lower_bound, validated_upper_bound = validate_increasing_parameter_bounds(
        parameter_name="search parameter",
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        owner_name="bounded likelihood search",
    )
    if tolerance <= 0.0:
        raise ValueError("bounded likelihood search tolerance must be positive")
    if max_iterations < 1:
        raise ValueError("bounded likelihood search iterations must be positive")

    phi = (math.sqrt(5.0) - 1.0) / 2.0
    cache: dict[float, tuple[_PayloadT, float]] = {}

    def evaluate_cached(parameter_value: float) -> tuple[_PayloadT, float]:
        validated_parameter_value = validate_parameter_within_bounds(
            parameter_name="search parameter",
            value=parameter_value,
            lower_bound=validated_lower_bound,
            upper_bound=validated_upper_bound,
            owner_name="bounded likelihood search",
        )
        cache_key = round(validated_parameter_value, 12)
        cached = cache.get(cache_key)
        if cached is not None:
            return cached
        candidate = evaluate(validated_parameter_value)
        cache[cache_key] = candidate
        return candidate

    left = validated_upper_bound - (
        phi * (validated_upper_bound - validated_lower_bound)
    )
    right = validated_lower_bound + (
        phi * (validated_upper_bound - validated_lower_bound)
    )
    left_payload, left_objective = evaluate_cached(left)
    right_payload, right_objective = evaluate_cached(right)
    iteration = 0

    while (
        validated_upper_bound - validated_lower_bound
    ) > tolerance and iteration < max_iterations:
        if left_objective < right_objective:
            validated_lower_bound = left
            left = right
            left_payload = right_payload
            left_objective = right_objective
            right = validated_lower_bound + (
                phi * (validated_upper_bound - validated_lower_bound)
            )
            right_payload, right_objective = evaluate_cached(right)
        else:
            validated_upper_bound = right
            right = left
            right_payload = left_payload
            right_objective = left_objective
            left = validated_upper_bound - (
                phi * (validated_upper_bound - validated_lower_bound)
            )
            left_payload, left_objective = evaluate_cached(left)
        iteration += 1

    midpoint = (validated_lower_bound + validated_upper_bound) / 2.0
    midpoint_payload, midpoint_objective = evaluate_cached(midpoint)
    best_parameter, best_payload, best_objective = sorted(
        [
            (left, left_payload, left_objective),
            (right, right_payload, right_objective),
            (midpoint, midpoint_payload, midpoint_objective),
        ],
        key=lambda item: (-item[2], item[0]),
    )[0]
    return BoundedLikelihoodSearchResult(
        parameter_value=best_parameter,
        payload=best_payload,
        objective_value=best_objective,
        function_evaluation_count=len(cache),
        converged=(validated_upper_bound - validated_lower_bound) <= tolerance,
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
            function_evaluation_count += search_result.function_evaluation_count
            if search_result.objective_value > (
                current_objective + improvement_tolerance
            ):
                current_values[name] = search_result.parameter_value
                current_payload = search_result.payload
                current_objective = search_result.objective_value
                improved = True
        if not improved:
            converged = True
            break

    return BoundedCoordinateLikelihoodSearchResult(
        parameter_values=dict(current_values),
        payload=current_payload,
        objective_value=current_objective,
        function_evaluation_count=function_evaluation_count,
        optimization_pass_count=optimization_pass_count,
        converged=converged,
    )
