from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from statistics import mean


@dataclass(slots=True)
class TraceSeriesSummary:
    path: Path
    parameter: str
    sample_count: int
    effective_sample_size: float
    mean: float
    minimum: float
    maximum: float
    first_half_mean: float
    second_half_mean: float
    standardized_mean_shift: float


@dataclass(slots=True)
class TraceConvergenceWarning:
    path: Path
    parameter: str
    code: str
    message: str
    observed_value: float
    threshold: float


@dataclass(slots=True)
class TraceConvergenceReport:
    path: Path
    sample_count: int
    ess_threshold: float
    mean_shift_threshold: float
    converged: bool
    series: list[TraceSeriesSummary]
    warnings: list[TraceConvergenceWarning]


def autocorrelation(series: list[float], lag: int) -> float:
    """Compute lag-k autocorrelation for one numeric trace series."""
    sample_count = len(series)
    average = mean(series)
    denominator = sum((value - average) ** 2 for value in series)
    if denominator == 0.0:
        return 0.0
    numerator = sum((series[index] - average) * (series[index + lag] - average) for index in range(sample_count - lag))
    return numerator / denominator


def effective_sample_size(series: list[float]) -> float:
    """Estimate effective sample size from positive autocorrelation decay."""
    sample_count = len(series)
    if sample_count < 3:
        return float(sample_count)
    rho_sum = 0.0
    for lag in range(1, sample_count):
        rho = autocorrelation(series, lag)
        if rho <= 0:
            break
        rho_sum += rho
    return round(sample_count / (1.0 + 2.0 * rho_sum), 6)


def standardized_mean_shift(series: list[float]) -> float:
    """Measure first-half versus second-half drift scaled by the whole-series spread."""
    sample_count = len(series)
    if sample_count < 4:
        return 0.0
    midpoint = sample_count // 2
    left = series[:midpoint]
    right = series[midpoint:]
    if not left or not right:
        return 0.0
    average = mean(series)
    variance = sum((value - average) ** 2 for value in series) / sample_count
    if variance <= 0.0:
        return 0.0
    return round(abs(mean(left) - mean(right)) / (variance ** 0.5), 6)


def summarize_trace_convergence(
    *,
    path: Path,
    rows: list[dict[str, float]],
    columns: list[str],
    ess_threshold: float = 200.0,
    mean_shift_threshold: float = 0.5,
) -> TraceConvergenceReport:
    """Summarize convergence risks for a numeric MCMC trace table."""
    if not rows:
        raise ValueError("trace convergence requires at least one sampled row")
    series_summaries: list[TraceSeriesSummary] = []
    warnings: list[TraceConvergenceWarning] = []
    for parameter in columns:
        values = [row[parameter] for row in rows if parameter in row]
        if not values:
            continue
        midpoint = len(values) // 2
        left = values[:midpoint] or values
        right = values[midpoint:] or values
        ess = effective_sample_size(values)
        mean_shift = standardized_mean_shift(values)
        series_summaries.append(
            TraceSeriesSummary(
                path=path,
                parameter=parameter,
                sample_count=len(values),
                effective_sample_size=ess,
                mean=round(mean(values), 6),
                minimum=min(values),
                maximum=max(values),
                first_half_mean=round(mean(left), 6),
                second_half_mean=round(mean(right), 6),
                standardized_mean_shift=mean_shift,
            )
        )
        if ess < ess_threshold:
            warnings.append(
                TraceConvergenceWarning(
                    path=path,
                    parameter=parameter,
                    code="low-ess",
                    message="effective sample size is below the requested threshold",
                    observed_value=ess,
                    threshold=ess_threshold,
                )
            )
        if mean_shift > mean_shift_threshold:
            warnings.append(
                TraceConvergenceWarning(
                    path=path,
                    parameter=parameter,
                    code="mean-drift",
                    message="first-half and second-half trace means differ more than the allowed threshold",
                    observed_value=mean_shift,
                    threshold=mean_shift_threshold,
                )
            )
    return TraceConvergenceReport(
        path=path,
        sample_count=len(rows),
        ess_threshold=ess_threshold,
        mean_shift_threshold=mean_shift_threshold,
        converged=not warnings,
        series=series_summaries,
        warnings=warnings,
    )
