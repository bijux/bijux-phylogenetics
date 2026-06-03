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


def summarize_metropolis_hastings_trace_posterior_intervals(
    *,
    chain_report: MetropolisHastingsRunReport,
    mass_fraction: float = 0.95,
) -> MetropolisHastingsTracePosteriorIntervalReport:
    """Summarize HPD and equal-tail intervals per scalar parameter for one chain."""
    if not isinstance(chain_report, MetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "metropolis-hastings trace posterior interval summary requires one MetropolisHastingsRunReport",
            code="trace_posterior_interval_chain_report_type_invalid",
        )
    validated_mass_fraction = _validate_mass_fraction(mass_fraction)
    parameter_names = _resolve_scalar_parameter_names(chain_report)
    parameter_rows = []
    for parameter_name in parameter_names:
        values = _extract_scalar_parameter_trace(
            chain_report=chain_report,
            parameter_name=parameter_name,
        )
        hpd_interval = compute_highest_posterior_density_interval(
            values,
            mass_fraction=validated_mass_fraction,
        )
        equal_tail_interval = compute_equal_tail_interval(
            values,
            mass_fraction=validated_mass_fraction,
        )
        parameter_rows.append(
            TracePosteriorIntervalRow(
                parameter_name=parameter_name,
                sample_count=len(values),
                mass_fraction=validated_mass_fraction,
                hpd_lower_bound=hpd_interval.lower_bound,
                hpd_upper_bound=hpd_interval.upper_bound,
                equal_tail_lower_bound=equal_tail_interval[0],
                equal_tail_upper_bound=equal_tail_interval[1],
            )
        )
    return MetropolisHastingsTracePosteriorIntervalReport(
        sample_every=chain_report.sample_every,
        parameter_rows=parameter_rows,
    )


def summarize_independent_metropolis_hastings_trace_posterior_intervals(
    *,
    run_report: IndependentMetropolisHastingsRunReport,
    mass_fraction: float = 0.95,
) -> IndependentMetropolisHastingsTracePosteriorIntervalReport:
    """Summarize posterior intervals across named independent chains."""
    if not isinstance(run_report, IndependentMetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "independent metropolis-hastings trace posterior interval summary requires one IndependentMetropolisHastingsRunReport",
            code="trace_posterior_interval_independent_run_report_type_invalid",
        )
    return IndependentMetropolisHastingsTracePosteriorIntervalReport(
        chain_reports=[
            IndependentMetropolisHastingsChainTracePosteriorIntervalReport(
                chain_name=chain_report.chain_name,
                posterior_interval_report=(
                    summarize_metropolis_hastings_trace_posterior_intervals(
                        chain_report=chain_report.chain_report,
                        mass_fraction=mass_fraction,
                    )
                ),
            )
            for chain_report in run_report.chain_reports
        ]
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


def _resolve_scalar_parameter_names(
    chain_report: MetropolisHastingsRunReport,
) -> list[str]:
    sampled_states = list(chain_report.sampled_states)
    if not sampled_states:
        raise PhylogeneticsError(
            "metropolis-hastings trace posterior interval summary requires at least one sampled state",
            code="trace_posterior_interval_sampled_states_empty",
        )
    first_parameter_names = set(sampled_states[0].model_parameters.scalar_parameters)
    if not first_parameter_names:
        raise PhylogeneticsError(
            "metropolis-hastings trace posterior interval summary requires at least one scalar model parameter",
            code="trace_posterior_interval_scalar_parameters_empty",
        )
    for sampled_state in sampled_states[1:]:
        parameter_names = set(sampled_state.model_parameters.scalar_parameters)
        if parameter_names != first_parameter_names:
            raise PhylogeneticsError(
                "metropolis-hastings trace posterior interval summary requires the same scalar parameter set across sampled states",
                code="trace_posterior_interval_scalar_parameter_set_inconsistent",
                details={
                    "expected_scalar_parameters": sorted(first_parameter_names),
                    "observed_scalar_parameters": sorted(parameter_names),
                },
            )
    return sorted(first_parameter_names)


def _extract_scalar_parameter_trace(
    *,
    chain_report: MetropolisHastingsRunReport,
    parameter_name: str,
) -> list[float]:
    return [
        float(sampled_state.model_parameters.scalar_parameters[parameter_name])
        for sampled_state in chain_report.sampled_states
    ]


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


def _compute_linear_sample_quantile(
    values: Sequence[float], probability: float
) -> float:
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
    return (ordered_values[lower_index] * lower_weight) + (
        ordered_values[upper_index] * upper_weight
    )
