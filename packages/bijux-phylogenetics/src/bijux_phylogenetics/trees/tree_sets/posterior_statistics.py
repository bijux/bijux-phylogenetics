from __future__ import annotations

import math


def effective_sample_size(values: list[float]) -> float | None:
    count = len(values)
    if count == 0:
        return None
    if count == 1:
        return 1.0
    mean = sum(values) / count
    centered = [value - mean for value in values]
    variance = sum(value * value for value in centered) / count
    if variance <= 0.0:
        return float(count)
    autocorrelation_sum = 0.0
    for lag in range(1, count):
        numerator = sum(
            centered[index] * centered[index + lag] for index in range(count - lag)
        )
        denominator = (count - lag) * variance
        rho = numerator / denominator if denominator > 0.0 else 0.0
        if rho <= 0.0:
            break
        autocorrelation_sum += rho
    tau = 1.0 + (2.0 * autocorrelation_sum)
    return round(max(1.0, min(float(count), count / tau)), 15)


def highest_posterior_density_interval(
    values: list[float],
    *,
    mass: float,
) -> tuple[float, float]:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0], ordered[0]
    window_width = max(1, math.ceil(mass * len(ordered)))
    best_start = 0
    best_width = ordered[window_width - 1] - ordered[0]
    for start in range(1, len(ordered) - window_width + 1):
        width = ordered[start + window_width - 1] - ordered[start]
        if width < best_width:
            best_width = width
            best_start = start
    return ordered[best_start], ordered[best_start + window_width - 1]
