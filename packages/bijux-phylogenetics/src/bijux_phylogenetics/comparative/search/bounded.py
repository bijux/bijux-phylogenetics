from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math
from typing import Generic, TypeVar

from bijux_phylogenetics.runtime.errors import ComparativeMethodError

_PayloadT = TypeVar("_PayloadT")


@dataclass(frozen=True, slots=True)
class BoundedSearchControls:
    """Deterministic controls for one bounded one-parameter likelihood search."""

    coarse_grid_point_count: int = 81
    fine_grid_point_count: int = 81
    initial_parameter_value: float | None = None
    refinement_start_count: int = 3


@dataclass(slots=True)
class BoundedSearchProfileRow:
    """One evaluated point on a bounded objective surface."""

    parameter_value: float
    objective_value: float


@dataclass(slots=True)
class BoundedSearchDiagnostics:
    """Shared audit surface for bounded one-parameter maximization."""

    optimizer_name: str
    parameter_search_strategy: str
    lower_bound: float
    upper_bound: float
    starting_parameter_policy: str
    starting_parameter_value: float
    starting_parameter_objective_value: float
    coarse_grid_point_count: int
    fine_grid_point_count: int
    refinement_start_count: int
    function_evaluation_count: int
    coarse_best_parameter: float
    coarse_best_objective_value: float
    fine_search_start: float
    fine_search_stop: float
    converged: bool
    hit_lower_boundary: bool
    hit_upper_boundary: bool


@dataclass(slots=True)
class BoundedSearchResult(Generic[_PayloadT]):
    """Best payload and diagnostics from one bounded one-parameter search."""

    parameter_value: float
    payload: _PayloadT
    objective_value: float
    diagnostics: BoundedSearchDiagnostics
    profile_rows: list[BoundedSearchProfileRow]


def run_bounded_maximization(
    *,
    lower_bound: float,
    upper_bound: float,
    controls: BoundedSearchControls,
    evaluate: Callable[[float], tuple[_PayloadT, float]],
    optimizer_name: str = "governed-bounded-multi-start-search",
    parameter_search_strategy: str = "bounded-coarse-grid-with-multi-start-golden-section-refinement",
) -> BoundedSearchResult[_PayloadT]:
    """Maximize one deterministic bounded objective with governed multi-start refinement."""
    if upper_bound <= lower_bound:
        raise ComparativeMethodError("parameter bounds must be strictly increasing")
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

    coarse_grid = _linspace(lower_bound, upper_bound, controls.coarse_grid_point_count)
    coarse_candidates, starting_parameter_policy, starting_parameter_value = (
        _ordered_coarse_candidates(
            coarse_grid,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            initial_parameter_value=controls.initial_parameter_value,
        )
    )
    cache: dict[float, tuple[_PayloadT, float]] = {}

    def evaluate_cached(parameter_value: float) -> tuple[_PayloadT, float]:
        cache_key = round(parameter_value, 12)
        cached = cache.get(cache_key)
        if cached is not None:
            return cached
        candidate = evaluate(parameter_value)
        cache[cache_key] = candidate
        return candidate

    profile_by_parameter: dict[float, float] = {}
    coarse_rows: list[tuple[float, float]] = []
    coarse_results: list[tuple[float, _PayloadT, float]] = []
    best_parameter = coarse_candidates[0]
    best_payload, best_objective_value = evaluate_cached(best_parameter)
    starting_parameter_objective_value = best_objective_value
    profile_by_parameter[round(best_parameter, 12)] = best_objective_value
    coarse_rows.append((best_parameter, best_objective_value))
    coarse_results.append((best_parameter, best_payload, best_objective_value))
    coarse_best_parameter = best_parameter
    coarse_best_objective_value = best_objective_value
    for candidate in coarse_candidates[1:]:
        payload, objective_value = evaluate_cached(candidate)
        profile_by_parameter[round(candidate, 12)] = objective_value
        coarse_rows.append((candidate, objective_value))
        coarse_results.append((candidate, payload, objective_value))
        if objective_value > best_objective_value:
            best_parameter = candidate
            best_payload = payload
            best_objective_value = objective_value
            coarse_best_parameter = candidate
            coarse_best_objective_value = objective_value

    sorted_coarse_grid = sorted(
        {round(value, 12): value for value in coarse_grid}.values()
    )
    winning_bracket = _search_bracket(
        coarse_grid=sorted_coarse_grid,
        anchor_parameter=best_parameter,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
    )
    ranked_coarse_results = sorted(
        coarse_results,
        key=lambda item: (-item[2], item[0]),
    )
    start_count = min(controls.refinement_start_count, len(ranked_coarse_results))
    for candidate_parameter, _, _ in ranked_coarse_results[:start_count]:
        bracket = _search_bracket(
            coarse_grid=sorted_coarse_grid,
            anchor_parameter=candidate_parameter,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
        )
        candidate_best_parameter, candidate_best_payload, candidate_best_objective = (
            _golden_section_refinement(
                lower_bound=bracket[0],
                upper_bound=bracket[1],
                fine_grid_point_count=controls.fine_grid_point_count,
                evaluate=evaluate_cached,
                profile_by_parameter=profile_by_parameter,
            )
        )
        if candidate_best_objective > best_objective_value:
            best_parameter = candidate_best_parameter
            best_payload = candidate_best_payload
            best_objective_value = candidate_best_objective
            winning_bracket = bracket

    fine_step = 0.0
    if controls.fine_grid_point_count > 1:
        fine_step = (winning_bracket[1] - winning_bracket[0]) / float(
            controls.fine_grid_point_count - 1
        )
    tolerance = max(abs(fine_step) / 2.0, 1e-9)
    return BoundedSearchResult(
        parameter_value=best_parameter,
        payload=best_payload,
        objective_value=best_objective_value,
        diagnostics=BoundedSearchDiagnostics(
            optimizer_name=optimizer_name,
            parameter_search_strategy=parameter_search_strategy,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            starting_parameter_policy=starting_parameter_policy,
            starting_parameter_value=starting_parameter_value,
            starting_parameter_objective_value=starting_parameter_objective_value,
            coarse_grid_point_count=len(coarse_grid),
            fine_grid_point_count=controls.fine_grid_point_count,
            refinement_start_count=start_count,
            function_evaluation_count=len(cache),
            coarse_best_parameter=coarse_best_parameter,
            coarse_best_objective_value=coarse_best_objective_value,
            fine_search_start=winning_bracket[0],
            fine_search_stop=winning_bracket[1],
            converged=True,
            hit_lower_boundary=math.isclose(
                best_parameter,
                lower_bound,
                rel_tol=0.0,
                abs_tol=tolerance,
            ),
            hit_upper_boundary=math.isclose(
                best_parameter,
                upper_bound,
                rel_tol=0.0,
                abs_tol=tolerance,
            ),
        ),
        profile_rows=[
            BoundedSearchProfileRow(
                parameter_value=parameter_value,
                objective_value=objective_value,
            )
            for parameter_value, objective_value in sorted(
                (
                    (parameter_value, objective_value)
                    for parameter_value, objective_value in profile_by_parameter.items()
                ),
                key=lambda item: item[0],
            )
        ],
    )


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
    """Maximize one bounded objective with one deterministic golden-section search."""
    if upper_bound <= lower_bound:
        raise ComparativeMethodError("parameter bounds must be strictly increasing")
    if tolerance <= 0.0:
        raise ComparativeMethodError("tolerance must be positive for bounded search")
    if max_iterations < 1:
        raise ComparativeMethodError(
            "max_iterations must be at least 1 for bounded search"
        )

    phi = (math.sqrt(5.0) - 1.0) / 2.0
    original_lower_bound = lower_bound
    original_upper_bound = upper_bound
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
    left_payload, left_objective_value = evaluate_cached(left)
    right_payload, right_objective_value = evaluate_cached(right)
    starting_parameter_value = left
    starting_parameter_objective_value = left_objective_value
    converged = False

    for _ in range(max_iterations):
        if abs(upper_bound - lower_bound) < tolerance:
            converged = True
            break
        if left_objective_value > right_objective_value:
            upper_bound = right
            right = left
            right_payload = left_payload
            right_objective_value = left_objective_value
            left = upper_bound - (phi * (upper_bound - lower_bound))
            left_payload, left_objective_value = evaluate_cached(left)
        else:
            lower_bound = left
            left = right
            left_payload = right_payload
            left_objective_value = right_objective_value
            right = lower_bound + (phi * (upper_bound - lower_bound))
            right_payload, right_objective_value = evaluate_cached(right)

    if left_objective_value >= right_objective_value:
        best_parameter = left
        best_payload = left_payload
        best_objective_value = left_objective_value
    else:
        best_parameter = right
        best_payload = right_payload
        best_objective_value = right_objective_value
    boundary_tolerance = tolerance
    return BoundedSearchResult(
        parameter_value=best_parameter,
        payload=best_payload,
        objective_value=best_objective_value,
        diagnostics=BoundedSearchDiagnostics(
            optimizer_name=optimizer_name,
            parameter_search_strategy=parameter_search_strategy,
            lower_bound=original_lower_bound,
            upper_bound=original_upper_bound,
            starting_parameter_policy="left-golden-interior-first-evaluation",
            starting_parameter_value=starting_parameter_value,
            starting_parameter_objective_value=starting_parameter_objective_value,
            coarse_grid_point_count=2,
            fine_grid_point_count=0,
            refinement_start_count=1,
            function_evaluation_count=len(cache),
            coarse_best_parameter=best_parameter,
            coarse_best_objective_value=best_objective_value,
            fine_search_start=lower_bound,
            fine_search_stop=upper_bound,
            converged=converged,
            hit_lower_boundary=math.isclose(
                best_parameter,
                original_lower_bound,
                rel_tol=0.0,
                abs_tol=boundary_tolerance,
            ),
            hit_upper_boundary=math.isclose(
                best_parameter,
                original_upper_bound,
                rel_tol=0.0,
                abs_tol=boundary_tolerance,
            ),
        ),
        profile_rows=[
            BoundedSearchProfileRow(
                parameter_value=parameter_value,
                objective_value=objective_value,
            )
            for parameter_value, (_, objective_value) in sorted(
                cache.items(),
                key=lambda item: item[0],
            )
        ],
    )


def _ordered_coarse_candidates(
    coarse_candidates: list[float],
    *,
    lower_bound: float,
    upper_bound: float,
    initial_parameter_value: float | None,
) -> tuple[list[float], str, float]:
    if initial_parameter_value is None:
        return coarse_candidates, "lower-bound-first-evaluation", coarse_candidates[0]
    if initial_parameter_value < lower_bound or initial_parameter_value > upper_bound:
        raise ComparativeMethodError(
            "initial_parameter_value must fall within the declared bounded search interval"
        )
    ordered = [initial_parameter_value]
    ordered.extend(
        candidate
        for candidate in coarse_candidates
        if not math.isclose(
            candidate,
            initial_parameter_value,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    )
    return ordered, "user-provided-first-evaluation", initial_parameter_value


def _linspace(start: float, stop: float, count: int) -> list[float]:
    if count < 2:
        return [start]
    step = (stop - start) / float(count - 1)
    return [start + (step * index) for index in range(count)]


def _search_bracket(
    *,
    coarse_grid: list[float],
    anchor_parameter: float,
    lower_bound: float,
    upper_bound: float,
) -> tuple[float, float]:
    best_index = min(
        range(len(coarse_grid)),
        key=lambda index: abs(coarse_grid[index] - anchor_parameter),
    )
    left = coarse_grid[max(0, best_index - 1)]
    right = coarse_grid[min(len(coarse_grid) - 1, best_index + 1)]
    return max(left, lower_bound), min(right, upper_bound)


def _golden_section_refinement(
    *,
    lower_bound: float,
    upper_bound: float,
    fine_grid_point_count: int,
    evaluate: Callable[[float], tuple[_PayloadT, float]],
    profile_by_parameter: dict[float, float],
) -> tuple[float, _PayloadT, float]:
    if math.isclose(lower_bound, upper_bound, rel_tol=0.0, abs_tol=1e-15):
        payload, objective_value = evaluate(lower_bound)
        profile_by_parameter[round(lower_bound, 12)] = objective_value
        return lower_bound, payload, objective_value
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    left = upper_bound - phi * (upper_bound - lower_bound)
    right = lower_bound + phi * (upper_bound - lower_bound)
    left_payload, left_objective = evaluate(left)
    right_payload, right_objective = evaluate(right)
    profile_by_parameter[round(left, 12)] = left_objective
    profile_by_parameter[round(right, 12)] = right_objective
    for _ in range(max(fine_grid_point_count - 2, 0)):
        if abs(upper_bound - lower_bound) < 1e-6:
            break
        if left_objective > right_objective:
            upper_bound = right
            right = left
            right_payload = left_payload
            right_objective = left_objective
            left = upper_bound - phi * (upper_bound - lower_bound)
            left_payload, left_objective = evaluate(left)
            profile_by_parameter[round(left, 12)] = left_objective
            continue
        lower_bound = left
        left = right
        left_payload = right_payload
        left_objective = right_objective
        right = lower_bound + phi * (upper_bound - lower_bound)
        right_payload, right_objective = evaluate(right)
        profile_by_parameter[round(right, 12)] = right_objective
    if left_objective >= right_objective:
        return left, left_payload, left_objective
    return right, right_payload, right_objective
