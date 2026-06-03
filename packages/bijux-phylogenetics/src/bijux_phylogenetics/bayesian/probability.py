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


def normalize_log_probabilities(log_values: Sequence[float]) -> tuple[float, ...]:
    """Shift log probabilities so their exponentiated values sum to one."""
    if not log_values:
        raise ValueError("log probability normalization requires at least one value")

    validated_values = tuple(_validate_log_probability(value) for value in log_values)
    total_log_probability = logsumexp(validated_values)
    if total_log_probability == -math.inf:
        raise ValueError(
            "log probability normalization requires at least one finite value"
        )
    return tuple(value - total_log_probability for value in validated_values)


def compare_log_probabilities(
    left: float,
    right: float,
    *,
    tolerance: float = 0.0,
) -> int:
    """Compare two log probabilities without exponentiating them."""
    left_value = _validate_log_probability(left)
    right_value = _validate_log_probability(right)
    if math.isnan(tolerance) or math.isinf(tolerance) or tolerance < 0.0:
        raise ValueError(
            "log probability comparison tolerance must be finite and non-negative"
        )
    if left_value == right_value:
        return 0
    if left_value > right_value + tolerance:
        return 1
    if right_value > left_value + tolerance:
        return -1
    return 0


def _validate_log_probability(value: float) -> float:
    if math.isnan(value):
        raise ValueError("log probabilities must not be NaN")
    if value == math.inf:
        raise ValueError("log probabilities must not be positive infinity")
    return value


__all__ = [
    "compare_log_probabilities",
    "log_probability_add",
    "logsumexp",
    "normalize_log_probabilities",
]
