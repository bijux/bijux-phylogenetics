from __future__ import annotations

import pytest

from bijux_phylogenetics.core.bounded_parameter_validation import (
    validate_increasing_parameter_bounds,
    validate_parameter_within_bounds,
    validate_positive_parameter_bounds,
)


def test_validate_increasing_parameter_bounds_returns_finite_interval() -> None:
    assert validate_increasing_parameter_bounds(
        parameter_name="alpha",
        lower_bound=0.05,
        upper_bound=20.0,
        owner_name="bounded parameter validation test",
    ) == (0.05, 20.0)


def test_validate_increasing_parameter_bounds_rejects_nonincreasing_interval() -> None:
    with pytest.raises(
        ValueError,
        match="bounded parameter validation test requires strictly increasing bounds",
    ):
        validate_increasing_parameter_bounds(
            parameter_name="alpha",
            lower_bound=1.0,
            upper_bound=1.0,
            owner_name="bounded parameter validation test",
        )


def test_validate_parameter_within_bounds_rejects_out_of_range_value() -> None:
    with pytest.raises(
        ValueError,
        match=(
            "bounded parameter validation test requires 'alpha' to lie within "
            r"\[0.05, 20.0\]"
        ),
    ):
        validate_parameter_within_bounds(
            parameter_name="alpha",
            value=21.0,
            lower_bound=0.05,
            upper_bound=20.0,
            owner_name="bounded parameter validation test",
        )


def test_validate_positive_parameter_bounds_rejects_nonpositive_lower_bound() -> None:
    with pytest.raises(
        ValueError,
        match=(
            "bounded parameter validation test requires a positive lower bound for "
            "'alpha'"
        ),
    ):
        validate_positive_parameter_bounds(
            parameter_name="alpha",
            lower_bound=0.0,
            upper_bound=20.0,
            owner_name="bounded parameter validation test",
        )
