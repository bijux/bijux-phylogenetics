from __future__ import annotations

from decimal import Decimal, localcontext
import math

import pytest

from bijux_phylogenetics.bayesian.probability import (
    log_probability_add,
    logsumexp,
)


def test_logsumexp_matches_high_precision_underflow_fixture() -> None:
    log_values = (-1000.0, -1001.0, -1002.0)

    result = logsumexp(log_values)

    with localcontext() as context:
        context.prec = 80
        expected = (
            Decimal("-1000")
            + (
                Decimal(1)
                + Decimal("-1").exp()
                + Decimal("-2").exp()
            ).ln()
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
