from __future__ import annotations

from collections.abc import Sequence
import math


def logsumexp(log_values: Sequence[float]) -> float:
    """Return log(sum(exp(log_values))) without underflow."""
    if not log_values:
        raise ValueError("logsumexp requires at least one log value")

    validated_values = tuple(_validate_log_probability(value) for value in log_values)
    finite_values = [value for value in validated_values if value != -math.inf]
    if not finite_values:
        return -math.inf

    offset = max(finite_values)
    scaled_total = sum(math.exp(value - offset) for value in finite_values)
    return offset + math.log(scaled_total)


def log_probability_add(left: float, right: float) -> float:
    """Return log(exp(left) + exp(right)) without leaving log space."""
    return logsumexp((left, right))


def _validate_log_probability(value: float) -> float:
    if math.isnan(value):
        raise ValueError("log probabilities must not be NaN")
    if value == math.inf:
        raise ValueError("log probabilities must not be positive infinity")
    return value


__all__ = [
    "log_probability_add",
    "logsumexp",
]
