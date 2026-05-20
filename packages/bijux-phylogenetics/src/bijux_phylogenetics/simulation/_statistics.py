from __future__ import annotations

from math import sqrt
from statistics import median


def _round_float(value: float) -> float:
    return round(float(value), 15)


def _mean(values: list[float]) -> float:
    return _round_float(sum(values) / len(values))


def _median(values: list[float]) -> float:
    return _round_float(float(median(values)))


def _population_standard_deviation(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean_value = sum(values) / len(values)
    variance = sum((value - mean_value) ** 2 for value in values) / len(values)
    return _round_float(variance**0.5)


def _sample_standard_deviation(values: list[float]) -> float:
    if len(values) <= 1:
        return 0.0
    center = _mean(values)
    return _round_float(
        sqrt(sum((value - center) ** 2 for value in values) / (len(values) - 1))
    )


def _sample_covariance(left: list[float], right: list[float]) -> float:
    if len(left) != len(right):
        raise ValueError("covariance inputs must have the same length")
    if len(left) <= 1:
        return 0.0
    left_center = _mean(left)
    right_center = _mean(right)
    return _round_float(
        sum(
            (left_value - left_center) * (right_value - right_center)
            for left_value, right_value in zip(left, right, strict=True)
        )
        / (len(left) - 1)
    )


def _sample_correlation(left: list[float], right: list[float]) -> float:
    left_standard_deviation = _sample_standard_deviation(left)
    right_standard_deviation = _sample_standard_deviation(right)
    if left_standard_deviation == 0.0 or right_standard_deviation == 0.0:
        return 0.0
    return _round_float(
        _sample_covariance(left, right)
        / (left_standard_deviation * right_standard_deviation)
    )
