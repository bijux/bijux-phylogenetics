from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import Generic, TypeVar

from .bounded import (
    BoundedSearchControls,
    BoundedSearchDiagnostics,
    run_bounded_maximization,
)

_PayloadT = TypeVar("_PayloadT")


@dataclass(slots=True)
class ComparativeBoundedScalarProfileRow:
    """One profile row from one comparative bounded scalar fit."""

    parameter_value: float
    objective_value: float
    delta_objective_value: float
    within_95_confidence_interval: bool


@dataclass(slots=True)
class ComparativeBoundedScalarFitResult(Generic[_PayloadT]):
    """Shared comparative bounded-scalar fit result with full optimizer diagnostics."""

    parameter_value: float
    payload: _PayloadT
    objective_value: float
    lower_95_confidence_interval: float | None
    upper_95_confidence_interval: float | None
    optimizer_diagnostics: BoundedSearchDiagnostics
    profile_rows: list[ComparativeBoundedScalarProfileRow]


def fit_bounded_scalar_parameter(
    *,
    lower_bound: float,
    upper_bound: float,
    controls: BoundedSearchControls,
    evaluate: Callable[[float], tuple[_PayloadT, float]],
    profile_parameter_values: Iterable[float],
    confidence_interval_drop: float,
    optimizer_name: str = "governed-two-stage-grid-search",
    parameter_search_strategy: str = "bounded-two-stage-grid-search",
) -> ComparativeBoundedScalarFitResult[_PayloadT]:
    """Fit one bounded scalar parameter and build its reviewer-facing likelihood profile."""
    evaluation_cache: dict[float, tuple[_PayloadT, float]] = {}

    def evaluate_cached(parameter_value: float) -> tuple[_PayloadT, float]:
        cache_key = round(parameter_value, 12)
        cached = evaluation_cache.get(cache_key)
        if cached is not None:
            return cached
        result = evaluate(parameter_value)
        evaluation_cache[cache_key] = result
        return result

    search_result = run_bounded_maximization(
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        controls=controls,
        evaluate=evaluate_cached,
        optimizer_name=optimizer_name,
        parameter_search_strategy=parameter_search_strategy,
    )
    profile_rows: list[ComparativeBoundedScalarProfileRow] = []
    best_objective_value = search_result.objective_value
    for parameter_value in profile_parameter_values:
        _, objective_value = evaluate_cached(parameter_value)
        delta_objective_value = max(0.0, best_objective_value - objective_value)
        profile_rows.append(
            ComparativeBoundedScalarProfileRow(
                parameter_value=parameter_value,
                objective_value=objective_value,
                delta_objective_value=delta_objective_value,
                within_95_confidence_interval=(
                    delta_objective_value <= confidence_interval_drop
                ),
            )
        )
    supported_parameter_values = [
        row.parameter_value for row in profile_rows if row.within_95_confidence_interval
    ]
    return ComparativeBoundedScalarFitResult(
        parameter_value=search_result.parameter_value,
        payload=search_result.payload,
        objective_value=search_result.objective_value,
        lower_95_confidence_interval=(
            None if not supported_parameter_values else min(supported_parameter_values)
        ),
        upper_95_confidence_interval=(
            None if not supported_parameter_values else max(supported_parameter_values)
        ),
        optimizer_diagnostics=search_result.diagnostics,
        profile_rows=profile_rows,
    )
