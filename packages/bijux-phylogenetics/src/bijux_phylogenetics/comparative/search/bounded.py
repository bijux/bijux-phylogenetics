from __future__ import annotations

from collections.abc import Callable
from typing import TypeVar

from bijux_phylogenetics.core.bounded_scalar_search import (
    BoundedScalarSearchControls,
    BoundedScalarSearchDiagnostics,
    BoundedScalarSearchProfileRow,
    BoundedScalarSearchResult,
    run_bounded_scalar_golden_section_maximization,
    run_bounded_scalar_maximization,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

_PayloadT = TypeVar("_PayloadT")

BoundedSearchControls = BoundedScalarSearchControls
BoundedSearchDiagnostics = BoundedScalarSearchDiagnostics
BoundedSearchProfileRow = BoundedScalarSearchProfileRow
BoundedSearchResult = BoundedScalarSearchResult


def run_bounded_maximization(
    *,
    lower_bound: float,
    upper_bound: float,
    controls: BoundedSearchControls,
    evaluate: Callable[[float], tuple[_PayloadT, float]],
    optimizer_name: str = "governed-bounded-multi-start-search",
    parameter_search_strategy: str = "bounded-coarse-grid-with-multi-start-golden-section-refinement",
) -> BoundedSearchResult[_PayloadT]:
    """Maximize one bounded objective with the shared scalar search surface."""
    try:
        return run_bounded_scalar_maximization(
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            controls=controls,
            evaluate=evaluate,
            optimizer_name=optimizer_name,
            parameter_search_strategy=parameter_search_strategy,
        )
    except ValueError as error:
        raise ComparativeMethodError(str(error)) from error


def run_bounded_golden_section_maximization(
    *,
    lower_bound: float,
    upper_bound: float,
    evaluate: Callable[[float], tuple[_PayloadT, float]],
    tolerance: float = 1e-9,
    max_iterations: int = 400,
    optimizer_name: str = "golden-section-search",
    parameter_search_strategy: str = "bounded-single-start-golden-section-search",
) -> BoundedSearchResult[_PayloadT]:
    """Maximize one bounded objective with the shared scalar golden-section search."""
    try:
        return run_bounded_scalar_golden_section_maximization(
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            evaluate=evaluate,
            tolerance=tolerance,
            max_iterations=max_iterations,
            optimizer_name=optimizer_name,
            parameter_search_strategy=parameter_search_strategy,
        )
    except ValueError as error:
        raise ComparativeMethodError(str(error)) from error


__all__ = [
    "BoundedSearchControls",
    "BoundedSearchDiagnostics",
    "BoundedSearchProfileRow",
    "BoundedSearchResult",
    "run_bounded_golden_section_maximization",
    "run_bounded_maximization",
]
