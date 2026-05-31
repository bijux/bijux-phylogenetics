from __future__ import annotations

from collections.abc import Iterable, Sequence
import math


def logsumexp(log_values: Iterable[float]) -> float:
    """Return log(sum(exp(log_values))) while preserving tiny likelihood terms."""
    values = [float(value) for value in log_values]
    if not values:
        raise ValueError("logsumexp requires at least one log value")
    maximum = max(values)
    if math.isinf(maximum) and maximum < 0.0:
        return float("-inf")
    if not math.isfinite(maximum):
        raise ValueError("logsumexp requires finite log values or -inf sentinels")
    total = 0.0
    for value in values:
        if math.isinf(value) and value < 0.0:
            continue
        if not math.isfinite(value):
            raise ValueError("logsumexp requires finite log values or -inf sentinels")
        total += math.exp(value - maximum)
    if total <= 0.0 or not math.isfinite(total):
        return float("-inf")
    return maximum + math.log(total)


def log_weighted_sum_exp(
    log_values: Sequence[float],
    *,
    weights: Sequence[float],
) -> float:
    """Return log(sum(weight_i * exp(log_value_i))) with zero-safe weighting."""
    if len(log_values) != len(weights):
        raise ValueError("log_weighted_sum_exp requires matching value and weight counts")
    weighted_logs: list[float] = []
    for log_value, weight in zip(log_values, weights, strict=True):
        normalized_weight = float(weight)
        if not math.isfinite(normalized_weight) or normalized_weight < 0.0:
            raise ValueError("log_weighted_sum_exp weights must be finite and nonnegative")
        if normalized_weight == 0.0:
            continue
        weighted_logs.append(math.log(normalized_weight) + float(log_value))
    if not weighted_logs:
        return float("-inf")
    return logsumexp(weighted_logs)
