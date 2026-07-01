from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from typing import Generic, TypeVar

from bijux_phylogenetics.core.bounded_parameter_optimization import (
    fit_profiled_bounded_parameter,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .bounded import (
    BoundedSearchControls,
    BoundedSearchDiagnostics,
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
    parameter_name: str = "search parameter",
    owner_name: str = "comparative bounded parameter fit",
    optimizer_name: str = "governed-two-stage-grid-search",
    parameter_search_strategy: str = "bounded-two-stage-grid-search",
) -> ComparativeBoundedScalarFitResult[_PayloadT]:
    """Fit one bounded scalar parameter and build its reviewer-facing likelihood profile."""
    try:
        fit = fit_profiled_bounded_parameter(
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            controls=controls,
            evaluate=evaluate,
            profile_parameter_values=list(profile_parameter_values),
            confidence_interval_drop=confidence_interval_drop,
            parameter_name=parameter_name,
            owner_name=owner_name,
            optimizer_name=optimizer_name,
            parameter_search_strategy=parameter_search_strategy,
        )
    except ValueError as error:
        raise ComparativeMethodError(str(error)) from error

    return ComparativeBoundedScalarFitResult(
        parameter_value=fit.parameter_value,
        payload=fit.payload,
        objective_value=fit.objective_value,
        lower_95_confidence_interval=fit.lower_confidence_interval,
        upper_95_confidence_interval=fit.upper_confidence_interval,
        optimizer_diagnostics=fit.optimizer_diagnostics,
        profile_rows=[
            ComparativeBoundedScalarProfileRow(
                parameter_value=row.parameter_value,
                objective_value=row.objective_value,
                delta_objective_value=row.delta_objective_value,
                within_95_confidence_interval=row.within_confidence_interval,
            )
            for row in fit.profile_rows
        ],
    )
