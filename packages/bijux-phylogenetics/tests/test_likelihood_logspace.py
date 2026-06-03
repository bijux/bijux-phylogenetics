from __future__ import annotations

from decimal import Decimal, getcontext
import math

from bijux_phylogenetics.phylo.likelihood.logspace import (
    log_weighted_sum_exp,
    logsumexp,
)


def test_logsumexp_matches_high_precision_reference_for_tiny_terms() -> None:
    log_values = [-1000.0, -1003.5, -1011.25]

    observed = logsumexp(log_values)
    expected = _high_precision_logsumexp(log_values)

    assert math.isfinite(observed)
    assert math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-12)


def test_log_weighted_sum_exp_ignores_zero_weighted_terms_and_order() -> None:
    log_values = [-980.0, -1000.0, -990.0]
    weights = [0.0, 0.25, 0.75]

    forward = log_weighted_sum_exp(log_values, weights=weights)
    reversed_observed = log_weighted_sum_exp(
        list(reversed(log_values)),
        weights=list(reversed(weights)),
    )
    expected = _high_precision_log_weighted_sum_exp(log_values, weights)

    assert math.isfinite(forward)
    assert math.isclose(forward, expected, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(reversed_observed, expected, rel_tol=0.0, abs_tol=1e-12)


def test_log_weighted_sum_exp_returns_negative_infinity_when_every_term_is_zero() -> (
    None
):
    assert log_weighted_sum_exp(
        [float("-inf"), float("-inf")],
        weights=[0.0, 0.0],
    ) == float("-inf")


def _high_precision_logsumexp(log_values: list[float]) -> float:
    return _high_precision_log_weighted_sum_exp(
        log_values,
        [1.0] * len(log_values),
    )


def _high_precision_log_weighted_sum_exp(
    log_values: list[float],
    weights: list[float],
) -> float:
    getcontext().prec = 90
    maximum = max(Decimal(str(value)) for value in log_values)
    total = Decimal("0")
    for log_value, weight in zip(log_values, weights, strict=True):
        decimal_weight = Decimal(str(weight))
        if decimal_weight == 0:
            continue
        total += decimal_weight * (Decimal(str(log_value)) - maximum).exp()
    if total == 0:
        return float("-inf")
    return float(maximum + total.ln())
