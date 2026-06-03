from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from bijux_phylogenetics.bayesian.independent_chains import (
    IndependentMetropolisHastingsRunReport,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import MetropolisHastingsRunReport
from bijux_phylogenetics.bayesian.trace_autocorrelation import (
    compute_trace_autocorrelation,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


@dataclass(frozen=True, slots=True)
class TraceEffectiveSampleSizeRow:
    """One autocorrelation-time ESS estimate for one scalar trace parameter."""

    parameter_name: str
    sample_count: int
    integrated_autocorrelation_time: float
    effective_sample_size: float
    last_positive_lag: int | None


@dataclass(frozen=True, slots=True)
class MetropolisHastingsTraceEffectiveSampleSizeReport:
    """One per-parameter ESS report for one native Metropolis-Hastings chain."""

    sample_every: int
    parameter_rows: list[TraceEffectiveSampleSizeRow]


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsChainTraceEffectiveSampleSizeReport:
    """One named chain ESS report inside one independent-chain collection."""

    chain_name: str
    effective_sample_size_report: MetropolisHastingsTraceEffectiveSampleSizeReport


@dataclass(frozen=True, slots=True)
class IndependentMetropolisHastingsTraceEffectiveSampleSizeReport:
    """One ESS report collection across named independent chains."""

    chain_reports: list[
        IndependentMetropolisHastingsChainTraceEffectiveSampleSizeReport
    ]


def compute_trace_integrated_autocorrelation_time(
    values: Sequence[float],
    *,
    maximum_lag: int | None = None,
) -> tuple[float, int | None]:
    """Compute integrated autocorrelation time from positive-lag decay."""
    validated_values = _validate_numeric_series(values)
    resolved_maximum_lag = _resolve_maximum_lag(
        sample_count=len(validated_values),
        maximum_lag=maximum_lag,
        owner_name="trace integrated autocorrelation time computation",
    )
    if resolved_maximum_lag == 0:
        return 1.0, None
    autocorrelation_sum = 0.0
    last_positive_lag: int | None = None
    for lag in range(1, resolved_maximum_lag + 1):
        rho = compute_trace_autocorrelation(validated_values, lag=lag)
        if rho <= 0.0:
            break
        autocorrelation_sum += rho
        last_positive_lag = lag
    return round(1.0 + (2.0 * autocorrelation_sum), 15), last_positive_lag


def compute_trace_effective_sample_size(
    values: Sequence[float],
    *,
    maximum_lag: int | None = None,
) -> float:
    """Compute effective sample size from integrated autocorrelation time."""
    validated_values = _validate_numeric_series(values)
    integrated_autocorrelation_time, _last_positive_lag = (
        compute_trace_integrated_autocorrelation_time(
            validated_values,
            maximum_lag=maximum_lag,
        )
    )
    return round(len(validated_values) / integrated_autocorrelation_time, 15)


def summarize_metropolis_hastings_trace_effective_sample_size(
    *,
    chain_report: MetropolisHastingsRunReport,
    maximum_lag: int | None = None,
) -> MetropolisHastingsTraceEffectiveSampleSizeReport:
    """Summarize autocorrelation-time ESS per scalar parameter for one native chain."""
    if not isinstance(chain_report, MetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "metropolis-hastings trace effective sample size summary requires one MetropolisHastingsRunReport",
            code="trace_effective_sample_size_chain_report_type_invalid",
        )
    parameter_names = _resolve_scalar_parameter_names(chain_report)
    parameter_rows = []
    for parameter_name in parameter_names:
        values = _extract_scalar_parameter_trace(
            chain_report=chain_report,
            parameter_name=parameter_name,
        )
        integrated_autocorrelation_time, last_positive_lag = (
            compute_trace_integrated_autocorrelation_time(
                values,
                maximum_lag=maximum_lag,
            )
        )
        parameter_rows.append(
            TraceEffectiveSampleSizeRow(
                parameter_name=parameter_name,
                sample_count=len(values),
                integrated_autocorrelation_time=integrated_autocorrelation_time,
                effective_sample_size=compute_trace_effective_sample_size(
                    values,
                    maximum_lag=maximum_lag,
                ),
                last_positive_lag=last_positive_lag,
            )
        )
    return MetropolisHastingsTraceEffectiveSampleSizeReport(
        sample_every=chain_report.sample_every,
        parameter_rows=parameter_rows,
    )


def summarize_independent_metropolis_hastings_trace_effective_sample_size(
    *,
    run_report: IndependentMetropolisHastingsRunReport,
    maximum_lag: int | None = None,
) -> IndependentMetropolisHastingsTraceEffectiveSampleSizeReport:
    """Summarize autocorrelation-time ESS across named independent chains."""
    if not isinstance(run_report, IndependentMetropolisHastingsRunReport):
        raise PhylogeneticsError(
            "independent metropolis-hastings trace effective sample size summary requires one IndependentMetropolisHastingsRunReport",
            code="trace_effective_sample_size_independent_run_report_type_invalid",
        )
    return IndependentMetropolisHastingsTraceEffectiveSampleSizeReport(
        chain_reports=[
            IndependentMetropolisHastingsChainTraceEffectiveSampleSizeReport(
                chain_name=chain_report.chain_name,
                effective_sample_size_report=(
                    summarize_metropolis_hastings_trace_effective_sample_size(
                        chain_report=chain_report.chain_report,
                        maximum_lag=maximum_lag,
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
            "trace effective sample size computation requires at least one sampled value",
            code="trace_effective_sample_size_series_empty",
        )
    for value in validated_values:
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise PhylogeneticsError(
                "trace effective sample size computation requires every sampled value to be numeric",
                code="trace_effective_sample_size_series_value_type_invalid",
            )
    return [float(value) for value in validated_values]


def _resolve_scalar_parameter_names(
    chain_report: MetropolisHastingsRunReport,
) -> list[str]:
    sampled_states = list(chain_report.sampled_states)
    if not sampled_states:
        raise PhylogeneticsError(
            "metropolis-hastings trace effective sample size summary requires at least one sampled state",
            code="trace_effective_sample_size_sampled_states_empty",
        )
    first_parameter_names = set(sampled_states[0].model_parameters.scalar_parameters)
    if not first_parameter_names:
        raise PhylogeneticsError(
            "metropolis-hastings trace effective sample size summary requires at least one scalar model parameter",
            code="trace_effective_sample_size_scalar_parameters_empty",
        )
    for sampled_state in sampled_states[1:]:
        parameter_names = set(sampled_state.model_parameters.scalar_parameters)
        if parameter_names != first_parameter_names:
            raise PhylogeneticsError(
                "metropolis-hastings trace effective sample size summary requires the same scalar parameter set across sampled states",
                code="trace_effective_sample_size_scalar_parameter_set_inconsistent",
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
    if sample_count == 1:
        if maximum_lag not in (None, 0):
            raise PhylogeneticsError(
                f"{owner_name} requires 'maximum_lag' to equal zero when only one sampled value is available",
                code="trace_effective_sample_size_maximum_lag_singleton_invalid",
                details={"maximum_lag": maximum_lag, "sample_count": sample_count},
            )
        return 0
    if maximum_lag is None:
        return sample_count - 1
    validated_maximum_lag = _validate_nonnegative_integer(
        value=maximum_lag,
        field_name="maximum_lag",
        owner_name=owner_name,
    )
    if validated_maximum_lag >= sample_count:
        raise PhylogeneticsError(
            f"{owner_name} requires 'maximum_lag' to be smaller than the trace length",
            code="trace_effective_sample_size_maximum_lag_out_of_range",
            details={
                "maximum_lag": validated_maximum_lag,
                "sample_count": sample_count,
            },
        )
    return validated_maximum_lag


def _validate_nonnegative_integer(
    *,
    value: int,
    field_name: str,
    owner_name: str,
) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be one integer",
            code="trace_effective_sample_size_integer_type_invalid",
            details={"field_name": field_name},
        )
    if value < 0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be nonnegative",
            code="trace_effective_sample_size_integer_negative",
            details={"field_name": field_name, "value": value},
        )
    return value
