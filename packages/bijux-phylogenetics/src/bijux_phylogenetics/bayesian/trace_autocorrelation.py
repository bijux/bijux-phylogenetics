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


def summarize_metropolis_hastings_trace_autocorrelation(
    *,
    chain_report: MetropolisHastingsRunReport,
    maximum_lag: int | None = None,
) -> MetropolisHastingsTraceAutocorrelationReport:
    """Summarize lag autocorrelation per scalar parameter for one native chain."""
    if not isinstance(chain_report, MetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "metropolis-hastings trace autocorrelation summary requires one MetropolisHastingsRunReport",
            code="trace_autocorrelation_chain_report_type_invalid",
        )
    scalar_parameter_names = _resolve_scalar_parameter_names(chain_report)
    sample_count = len(chain_report.sampled_states)
    validated_maximum_lag = _resolve_maximum_lag(
        sample_count=sample_count,
        maximum_lag=maximum_lag,
        owner_name="metropolis-hastings trace autocorrelation summary",
    )
    parameter_reports = [
        TraceAutocorrelationParameterReport(
            parameter_name=parameter_name,
            sample_count=sample_count,
            lag_rows=[
                TraceAutocorrelationLagRow(
                    parameter_name=parameter_name,
                    lag=lag,
                    autocorrelation=compute_trace_autocorrelation(
                        _extract_scalar_parameter_trace(
                            chain_report=chain_report,
                            parameter_name=parameter_name,
                        ),
                        lag=lag,
                    ),
                )
                for lag in range(1, validated_maximum_lag + 1)
            ],
        )
        for parameter_name in scalar_parameter_names
    ]
    return MetropolisHastingsTraceAutocorrelationReport(
        sample_every=chain_report.sample_every,
        parameter_reports=parameter_reports,
    )


def summarize_independent_metropolis_hastings_trace_autocorrelation(
    *,
    run_report: IndependentMetropolisHastingsRunReport,
    maximum_lag: int | None = None,
) -> IndependentMetropolisHastingsTraceAutocorrelationReport:
    """Summarize lag autocorrelation per scalar parameter across named chains."""
    if not isinstance(run_report, IndependentMetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "independent metropolis-hastings trace autocorrelation summary requires one IndependentMetropolisHastingsRunReport",
            code="trace_autocorrelation_independent_run_report_type_invalid",
        )
    return IndependentMetropolisHastingsTraceAutocorrelationReport(
        chain_reports=[
            IndependentMetropolisHastingsChainTraceAutocorrelationReport(
                chain_name=chain_report.chain_name,
                autocorrelation_report=summarize_metropolis_hastings_trace_autocorrelation(
                    chain_report=chain_report.chain_report,
                    maximum_lag=maximum_lag,
                ),
            )
            for chain_report in run_report.chain_reports
        ]
    )


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


def _resolve_scalar_parameter_names(
    chain_report: MetropolisHastingsRunReport,
) -> list[str]:
    sampled_states = list(chain_report.sampled_states)
    if len(sampled_states) < 2:
        raise PhylogeneticsError(
            "metropolis-hastings trace autocorrelation summary requires at least two sampled states",
            code="trace_autocorrelation_sampled_state_count_too_small",
            details={"sampled_state_count": len(sampled_states)},
        )
    first_parameter_names = set(sampled_states[0].model_parameters.scalar_parameters)
    if not first_parameter_names:
        raise PhylogeneticsError(
            "metropolis-hastings trace autocorrelation summary requires at least one scalar model parameter",
            code="trace_autocorrelation_scalar_parameters_empty",
        )
    for sampled_state in sampled_states[1:]:
        parameter_names = set(sampled_state.model_parameters.scalar_parameters)
        if parameter_names != first_parameter_names:
            raise PhylogeneticsError(
                "metropolis-hastings trace autocorrelation summary requires the same scalar parameter set across sampled states",
                code="trace_autocorrelation_scalar_parameter_set_inconsistent",
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


def _resolve_maximum_lag(
    *,
    sample_count: int,
    maximum_lag: int | None,
    owner_name: str,
) -> int:
    if maximum_lag is None:
        return sample_count - 1
    validated_maximum_lag = _validate_positive_integer(
        value=maximum_lag,
        field_name="maximum_lag",
        owner_name=owner_name,
    )
    if validated_maximum_lag >= sample_count:
        raise PhylogeneticsError(
            f"{owner_name} requires 'maximum_lag' to be smaller than the sampled trace length",
            code="trace_autocorrelation_maximum_lag_out_of_range",
            details={
                "maximum_lag": validated_maximum_lag,
                "sample_count": sample_count,
            },
        )
    return validated_maximum_lag


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
