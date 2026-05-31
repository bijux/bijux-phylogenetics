from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from statistics import mean

from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsRunReport,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import MetropolisHastingsRunReport
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


@dataclass(frozen=True, slots=True)
class TraceAutocorrelationLagRow:
    """One lag-specific autocorrelation estimate for one scalar trace parameter."""

    parameter_name: str
    lag: int
    autocorrelation: float


@dataclass(frozen=True, slots=True)
class TraceAutocorrelationParameterReport:
    """One autocorrelation report for one scalar trace parameter."""

    parameter_name: str
    sample_count: int
    lag_rows: list[TraceAutocorrelationLagRow]


@dataclass(frozen=True, slots=True)
class MetropolisHastingsTraceAutocorrelationReport:
    """One scalar-parameter autocorrelation report for one native chain."""

    sample_every: int
    parameter_reports: list[TraceAutocorrelationParameterReport]


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsChainTraceAutocorrelationReport:
    """One named chain autocorrelation report inside one multi-chain run."""

    chain_name: str
    autocorrelation_report: MetropolisHastingsTraceAutocorrelationReport


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsTraceAutocorrelationReport:
    """One chainwise scalar-parameter autocorrelation report over many chains."""

    chain_reports: list[IndependentMetropolisHastingsChainTraceAutocorrelationReport]


def compute_trace_autocorrelation(
    values: Sequence[float],
    *,
    lag: int,
) -> float:
    """Compute one lag autocorrelation for one scalar numeric trace."""
    validated_values = _validate_numeric_series(values)
    validated_lag = _validate_positive_integer(
        value=lag,
        field_name="lag",
        owner_name="trace autocorrelation computation",
    )
    if validated_lag >= len(validated_values):
        raise PhylogeneticsError(
            "trace autocorrelation computation requires 'lag' to be smaller than the trace length",
            code="trace_autocorrelation_lag_out_of_range",
            details={
                "lag": validated_lag,
                "sample_count": len(validated_values),
            },
        )
    average = mean(validated_values)
    denominator = sum((value - average) ** 2 for value in validated_values)
    if denominator == 0.0:
        return 0.0
    numerator = sum(
        (validated_values[index] - average)
        * (validated_values[index + validated_lag] - average)
        for index in range(len(validated_values) - validated_lag)
    )
    return round(numerator / denominator, 15)


def _validate_numeric_series(values: Sequence[float]) -> list[float]:
    validated_values = list(values)
    if len(validated_values) < 2:
        raise PhylogeneticsError(
            "trace autocorrelation computation requires at least two sampled values",
            code="trace_autocorrelation_series_too_short",
            details={"sample_count": len(validated_values)},
        )
    for value in validated_values:
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise PhylogeneticsError(
                "trace autocorrelation computation requires every sampled value to be numeric",
                code="trace_autocorrelation_series_value_type_invalid",
            )
    return [float(value) for value in validated_values]


def _validate_positive_integer(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one integer",
            code="trace_autocorrelation_integer_type_invalid",
            details={"field_name": field_name},
        )
    if value <= 0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be positive",
            code="trace_autocorrelation_integer_not_positive",
            details={"field_name": field_name, "value": value},
        )
    return value
