from __future__ import annotations

from decimal import Decimal, localcontext
import math

import pytest

from bijux_phylogenetics.bayesian.probability import (
    compare_log_probabilities,
    log_probability_add,
    logsumexp,
    normalize_log_probabilities,
)


def test_logsumexp_matches_high_precision_underflow_fixture() -> None:
    log_values = (-1000.0, -1001.0, -1002.0)

    result = logsumexp(log_values)

    with localcontext() as context:
        context.prec = 80
        expected = (
            Decimal("-1000")
            + (Decimal(1) + Decimal("-1").exp() + Decimal("-2").exp()).ln()
        )
    assert math.isclose(result, float(expected), rel_tol=0.0, abs_tol=1e-12)


def test_logsumexp_returns_negative_infinity_for_zero_mass_distribution() -> None:
    result = logsumexp((-math.inf, -math.inf, -math.inf))

    assert result == -math.inf


def test_log_probability_add_matches_logsumexp_for_two_values() -> None:
    left = -1000.0
    right = -1002.5

    assert math.isclose(
        log_probability_add(left, right),
        logsumexp((left, right)),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_normalize_log_probabilities_matches_high_precision_fixture() -> None:
    log_values = (-1000.0, -1001.0, -1002.0)

    normalized = normalize_log_probabilities(log_values)

    with localcontext() as context:
        context.prec = 80
        total = (
            Decimal("-1000").exp() + Decimal("-1001").exp() + Decimal("-1002").exp()
        ).ln()
        expected = tuple(float(Decimal(str(value)) - total) for value in log_values)
    assert all(
        math.isclose(observed, wanted, rel_tol=0.0, abs_tol=1e-12)
        for observed, wanted in zip(normalized, expected, strict=True)
    )
    assert math.isclose(
        sum(math.exp(value) for value in normalized),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_normalize_log_probabilities_preserves_zero_mass_entries() -> None:
    normalized = normalize_log_probabilities((-1000.0, -math.inf, -1003.0))

    assert normalized[1] == -math.inf
    assert math.isclose(
        sum(math.exp(value) for value in normalized if value != -math.inf),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_normalize_log_probabilities_rejects_all_zero_mass_inputs() -> None:
    with pytest.raises(ValueError, match="at least one finite value"):
        normalize_log_probabilities((-math.inf, -math.inf))


def test_compare_log_probabilities_orders_values_without_exponentiating() -> None:
    assert compare_log_probabilities(-1000.0, -1001.0) == 1
    assert compare_log_probabilities(-1001.0, -1000.0) == -1
    assert compare_log_probabilities(-math.inf, -math.inf) == 0


def test_compare_log_probabilities_respects_tolerance() -> None:
    assert compare_log_probabilities(-12.0, -12.05, tolerance=0.1) == 0
    assert compare_log_probabilities(-12.0, -12.2, tolerance=0.1) == 1


@pytest.mark.parametrize(
    ("log_values", "message"),
    [
        ((), "at least one log value"),
        ((math.nan,), "must not be NaN"),
        ((math.inf,), "must not be positive infinity"),
    ],
)
def test_logsumexp_rejects_invalid_inputs(
    log_values: tuple[float, ...],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        logsumexp(log_values)


@pytest.mark.parametrize("tolerance", (math.nan, math.inf, -0.1))
def test_compare_log_probabilities_rejects_invalid_tolerance(
    tolerance: float,
) -> None:
    with pytest.raises(
        ValueError,
        match="tolerance must be finite and non-negative",
    ):
        compare_log_probabilities(-1.0, -2.0, tolerance=tolerance)
