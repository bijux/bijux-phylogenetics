from __future__ import annotations

import math

from bijux_phylogenetics.phylo.likelihood.models import DiscreteGammaRateCategory


def build_discrete_gamma_rate_categories(
    *,
    alpha: float,
    category_count: int,
) -> list[DiscreteGammaRateCategory]:
    """Build equal-probability midpoint discrete-gamma rate categories with unit mean."""
    validated_alpha = validate_discrete_gamma_alpha(alpha)
    validated_category_count = validate_discrete_gamma_category_count(category_count)
    raw_rates = [
        _gamma_quantile(
            (category_index + 0.5) / validated_category_count,
            validated_alpha,
        )
        for category_index in range(validated_category_count)
    ]
    normalized_rates = _normalize_weighted_mean_to_one(raw_rates)
    weight = 1.0 / validated_category_count
    return [
        DiscreteGammaRateCategory(
            category_index=category_index + 1,
            rate=rate,
            weight=weight,
        )
        for category_index, rate in enumerate(normalized_rates)
    ]


def validate_discrete_gamma_alpha(alpha: float) -> float:
    if not math.isfinite(alpha) or alpha <= 0.0:
        raise ValueError("discrete-gamma alpha must be a finite positive value")
    return float(alpha)


def validate_discrete_gamma_category_count(category_count: int) -> int:
    if category_count < 1:
        raise ValueError("discrete-gamma category_count must be at least one")
    return int(category_count)


def _normalize_weighted_mean_to_one(raw_rates: list[float]) -> list[float]:
    mean_rate = sum(raw_rates) / len(raw_rates)
    if mean_rate <= 0.0 or not math.isfinite(mean_rate):
        raise ValueError("discrete-gamma midpoint rates must have positive finite mean")
    return [rate / mean_rate for rate in raw_rates]


def _gamma_quantile(probability: float, alpha: float) -> float:
    if probability <= 0.0:
        return 0.0
    if probability >= 1.0:
        return float("inf")
    lower = 0.0
    upper = 1.0
    while _gamma_cdf(upper, alpha) < probability:
        upper *= 2.0
        if upper > 1_000_000.0:
            raise ValueError("failed to bracket discrete-gamma quantile")
    for _ in range(80):
        midpoint = (lower + upper) / 2.0
        if _gamma_cdf(midpoint, alpha) < probability:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def _gamma_cdf(x_value: float, alpha: float) -> float:
    if x_value <= 0.0:
        return 0.0
    interval_count = max(1024, int(math.ceil(x_value * 256.0)))
    step = x_value / interval_count
    total = 0.0
    for interval_index in range(interval_count):
        sample = (interval_index + 0.5) * step
        total += _gamma_density(sample, alpha)
    cdf = step * total
    return min(max(cdf, 0.0), 1.0)


def _gamma_density(x_value: float, alpha: float) -> float:
    if x_value < 0.0:
        return 0.0
    scale_inverse = alpha
    log_density = (
        alpha * math.log(scale_inverse)
        - math.lgamma(alpha)
        + ((alpha - 1.0) * math.log(x_value))
        - (scale_inverse * x_value)
    )
    return math.exp(log_density)
