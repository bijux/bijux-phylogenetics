from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math
from typing import Generic, TypeVar

_PayloadT = TypeVar("_PayloadT")


@dataclass(slots=True)
class BoundedLikelihoodSearchResult(Generic[_PayloadT]):
    """Best result from one bounded one-parameter likelihood search."""

    parameter_value: float
    payload: _PayloadT
    objective_value: float
    function_evaluation_count: int
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
    if upper_bound <= lower_bound:
        raise ValueError("bounded likelihood search requires increasing parameter bounds")
    if tolerance <= 0.0:
        raise ValueError("bounded likelihood search tolerance must be positive")
    if max_iterations < 1:
        raise ValueError("bounded likelihood search iterations must be positive")

    phi = (math.sqrt(5.0) - 1.0) / 2.0
    cache: dict[float, tuple[_PayloadT, float]] = {}

    def evaluate_cached(parameter_value: float) -> tuple[_PayloadT, float]:
        cache_key = round(parameter_value, 12)
        cached = cache.get(cache_key)
        if cached is not None:
            return cached
        candidate = evaluate(parameter_value)
        cache[cache_key] = candidate
        return candidate

    left = upper_bound - (phi * (upper_bound - lower_bound))
    right = lower_bound + (phi * (upper_bound - lower_bound))
    left_payload, left_objective = evaluate_cached(left)
    right_payload, right_objective = evaluate_cached(right)
    iteration = 0

    while (upper_bound - lower_bound) > tolerance and iteration < max_iterations:
        if left_objective < right_objective:
            lower_bound = left
            left = right
            left_payload = right_payload
            left_objective = right_objective
            right = lower_bound + (phi * (upper_bound - lower_bound))
            right_payload, right_objective = evaluate_cached(right)
        else:
            upper_bound = right
            right = left
            right_payload = left_payload
            right_objective = left_objective
            left = upper_bound - (phi * (upper_bound - lower_bound))
            left_payload, left_objective = evaluate_cached(left)
        iteration += 1

    midpoint = (lower_bound + upper_bound) / 2.0
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
        converged=(upper_bound - lower_bound) <= tolerance,
    )
