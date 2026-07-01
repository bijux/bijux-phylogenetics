from __future__ import annotations

from collections.abc import Callable
from typing import TypeVar

from bijux_phylogenetics.core.bounded_parameter_optimization import (
    BoundedCoordinateOptimizationResult,
    BoundedParameterOptimizationResult,
    BoundedParameterSearchDiagnostics,
    run_bounded_coordinate_optimization,
    run_bounded_parameter_golden_section_maximization,
)

_PayloadT = TypeVar("_PayloadT")

BoundedLikelihoodSearchResult = BoundedParameterOptimizationResult
BoundedCoordinateLikelihoodSearchResult = BoundedCoordinateOptimizationResult
BoundedScalarSearchDiagnostics = BoundedParameterSearchDiagnostics


def run_bounded_likelihood_search(
    *,
    lower_bound: float,
    upper_bound: float,
    evaluate: Callable[[float], tuple[_PayloadT, float]],
    parameter_name: str = "search parameter",
    owner_name: str = "bounded likelihood search",
    tolerance: float = 1e-9,
    max_iterations: int = 400,
) -> BoundedLikelihoodSearchResult[_PayloadT]:
    """Maximize one bounded likelihood-like objective by shared scalar search."""
    return run_bounded_parameter_golden_section_maximization(
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        evaluate=evaluate,
        parameter_name=parameter_name,
        owner_name=owner_name,
        tolerance=tolerance,
        max_iterations=max_iterations,
    )


def run_bounded_coordinate_likelihood_search(
    *,
    initial_values: dict[str, float],
    bounds_by_name: dict[str, tuple[float, float]],
    evaluate: Callable[[dict[str, float]], tuple[_PayloadT, float]],
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
    owner_name: str = "coordinate likelihood search",
) -> BoundedCoordinateLikelihoodSearchResult[_PayloadT]:
    """Maximize one bounded likelihood-like objective by coordinate search."""
    return run_bounded_coordinate_optimization(
        initial_values=initial_values,
        bounds_by_name=bounds_by_name,
        evaluate=evaluate,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
        owner_name=owner_name,
    )
