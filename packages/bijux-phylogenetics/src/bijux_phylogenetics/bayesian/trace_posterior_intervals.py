from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsRunReport,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import MetropolisHastingsRunReport
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


@dataclass(frozen=True, slots=True)
class HighestPosteriorDensityInterval:
    """One shortest posterior interval spanning one requested posterior mass."""

    mass_fraction: float
    sample_count: int
    lower_bound: float
    upper_bound: float


@dataclass(frozen=True, slots=True)
class TracePosteriorIntervalRow:
    """One posterior-interval summary for one scalar trace parameter."""

    parameter_name: str
    sample_count: int
    mass_fraction: float
    hpd_lower_bound: float
    hpd_upper_bound: float
    equal_tail_lower_bound: float
    equal_tail_upper_bound: float


@dataclass(frozen=True, slots=True)
class MetropolisHastingsTracePosteriorIntervalReport:
    """One per-parameter posterior-interval report for one native chain."""

    sample_every: int
    parameter_rows: list[TracePosteriorIntervalRow]


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsChainTracePosteriorIntervalReport:
    """One named chain posterior-interval report inside one multi-chain run."""

    chain_name: str
    posterior_interval_report: MetropolisHastingsTracePosteriorIntervalReport


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsTracePosteriorIntervalReport:
    """One collection of posterior-interval reports across named chains."""

    chain_reports: list[IndependentMetropolisHastingsChainTracePosteriorIntervalReport]


def compute_highest_posterior_density_interval(
    values: Sequence[float],
    *,
    mass_fraction: float = 0.95,
) -> HighestPosteriorDensityInterval:
    """Compute the shortest sample interval spanning the requested posterior mass."""
    validated_values = _validate_numeric_series(values)
    validated_mass_fraction = _validate_mass_fraction(mass_fraction)
    ordered_values = sorted(validated_values)
    if len(ordered_values) == 1:
        return HighestPosteriorDensityInterval(
            mass_fraction=validated_mass_fraction,
            sample_count=1,
            lower_bound=ordered_values[0],
            upper_bound=ordered_values[0],
        )
    window_size = max(1, math.ceil(validated_mass_fraction * len(ordered_values)))
    if window_size >= len(ordered_values):
        return HighestPosteriorDensityInterval(
            mass_fraction=validated_mass_fraction,
            sample_count=len(ordered_values),
            lower_bound=ordered_values[0],
            upper_bound=ordered_values[-1],
        )
    best_start = 0
    best_width = ordered_values[window_size - 1] - ordered_values[0]
    for start in range(1, len(ordered_values) - window_size + 1):
        width = ordered_values[start + window_size - 1] - ordered_values[start]
        if width < best_width:
            best_start = start
            best_width = width
    return HighestPosteriorDensityInterval(
        mass_fraction=validated_mass_fraction,
        sample_count=len(ordered_values),
        lower_bound=ordered_values[best_start],
        upper_bound=ordered_values[best_start + window_size - 1],
    )


def compute_equal_tail_interval(
    values: Sequence[float],
    *,
    mass_fraction: float = 0.95,
) -> tuple[float, float]:
    """Compute the equal-tail interval spanning the requested posterior mass."""
    validated_values = _validate_numeric_series(values)
    validated_mass_fraction = _validate_mass_fraction(mass_fraction)
    if len(validated_values) == 1:
        return validated_values[0], validated_values[0]
    lower_tail_probability = (1.0 - validated_mass_fraction) / 2.0
    upper_tail_probability = 1.0 - lower_tail_probability
    return (
        _compute_linear_sample_quantile(validated_values, lower_tail_probability),
        _compute_linear_sample_quantile(validated_values, upper_tail_probability),
    )


def _validate_numeric_series(values: Sequence[float]) -> list[float]:
    validated_values = list(values)
    if not validated_values:
        raise PhylogeneticsError(
            "posterior interval computation requires at least one sampled value",
            code="trace_posterior_interval_series_empty",
        )
    for value in validated_values:
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise PhylogeneticsError(
                "posterior interval computation requires every sampled value to be numeric",
                code="trace_posterior_interval_series_value_type_invalid",
            )
    return [float(value) for value in validated_values]


def _validate_mass_fraction(value: float) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise PhylogeneticsError(
            "posterior interval computation requires 'mass_fraction' to be numeric",
            code="trace_posterior_interval_mass_fraction_type_invalid",
        )
    validated_value = float(value)
    if not 0.0 < validated_value <= 1.0:
        raise PhylogeneticsError(
            "posterior interval computation requires 'mass_fraction' to lie in the interval (0, 1]",
            code="trace_posterior_interval_mass_fraction_out_of_range",
            details={"mass_fraction": validated_value},
        )
    return validated_value


def _compute_linear_sample_quantile(values: Sequence[float], probability: float) -> float:
    ordered_values = sorted(values)
    if len(ordered_values) == 1:
        return ordered_values[0]
    fractional_index = (len(ordered_values) - 1) * probability
    lower_index = math.floor(fractional_index)
    upper_index = math.ceil(fractional_index)
    if lower_index == upper_index:
        return ordered_values[lower_index]
    upper_weight = fractional_index - lower_index
    lower_weight = 1.0 - upper_weight
    return (
        (ordered_values[lower_index] * lower_weight)
        + (ordered_values[upper_index] * upper_weight)
    )
