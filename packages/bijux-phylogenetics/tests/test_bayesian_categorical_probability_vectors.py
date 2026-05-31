from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.bayesian.probability_vectors import (
    build_categorical_probability_vector,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_build_categorical_probability_vector_preserves_valid_probabilities() -> None:
    vector = build_categorical_probability_vector(
        {"A": 0.7, "B": 0.2, "C": 0.1},
    )

    assert vector.states == ("A", "B", "C")
    assert vector.probabilities == (0.7, 0.2, 0.1)
    assert vector.total_probability == pytest.approx(1.0)
    assert vector.as_mapping() == {"A": 0.7, "B": 0.2, "C": 0.1}


def test_build_categorical_probability_vector_rejects_negative_probability() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        build_categorical_probability_vector({"A": 1.1, "B": -0.1})

    assert error_info.value.code == "categorical_probability_vector_value_negative"
    assert error_info.value.details["state"] == "B"
    assert error_info.value.details["probability"] == -0.1


@pytest.mark.parametrize("invalid_probability", (math.nan, math.inf, -math.inf))
def test_build_categorical_probability_vector_rejects_non_finite_probabilities(
    invalid_probability: float,
) -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        build_categorical_probability_vector({"A": invalid_probability, "B": 0.0})

    assert error_info.value.code == "categorical_probability_vector_value_not_finite"
    assert error_info.value.details["state"] == "A"
    if math.isnan(invalid_probability):
        assert math.isnan(error_info.value.details["probability"])
    else:
        assert error_info.value.details["probability"] == invalid_probability


def test_build_categorical_probability_vector_rejects_non_normalized_input() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        build_categorical_probability_vector({"A": 0.6, "B": 0.3})

    assert error_info.value.code == "categorical_probability_vector_not_normalized"
    assert error_info.value.details["total_probability"] == pytest.approx(0.9)
    assert error_info.value.details["expected_total"] == 1.0
    assert error_info.value.details["absolute_error"] == pytest.approx(0.1)


def test_build_categorical_probability_vector_fills_missing_expected_states() -> None:
    vector = build_categorical_probability_vector(
        {"A": 0.7, "B": 0.3},
        expected_states=("A", "B", "C"),
        missing_state_policy="fill-zero",
    )

    assert vector.states == ("A", "B", "C")
    assert vector.probabilities == (0.7, 0.3, 0.0)
    assert vector.missing_state_policy == "fill-zero"
    assert vector.probability_for("C") == 0.0


def test_build_categorical_probability_vector_rejects_missing_expected_states() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        build_categorical_probability_vector(
            {"A": 0.7, "B": 0.3},
            expected_states=("A", "B", "C"),
        )

    assert error_info.value.code == "categorical_probability_vector_missing_states"
    assert error_info.value.details["missing_states"] == ["C"]
    assert error_info.value.details["missing_state_policy"] == "reject"


def test_build_categorical_probability_vector_rejects_unexpected_states() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        build_categorical_probability_vector(
            {"A": 0.5, "B": 0.4, "X": 0.1},
            expected_states=("A", "B", "C"),
            missing_state_policy="fill-zero",
        )

    assert error_info.value.code == "categorical_probability_vector_unexpected_states"
    assert error_info.value.details["unexpected_states"] == ["X"]


def test_build_categorical_probability_vector_rejects_invalid_missing_state_policy() -> (
    None
):
    with pytest.raises(PhylogeneticsError) as error_info:
        build_categorical_probability_vector(
            {"A": 1.0},
            missing_state_policy="renormalize",
        )

    assert (
        error_info.value.code
        == "categorical_probability_vector_missing_state_policy_invalid"
    )
    assert error_info.value.details["missing_state_policy"] == "renormalize"


@pytest.mark.parametrize("normalization_tolerance", (math.nan, math.inf, -1e-9))
def test_build_categorical_probability_vector_rejects_invalid_tolerance(
    normalization_tolerance: float,
) -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        build_categorical_probability_vector(
            {"A": 1.0},
            normalization_tolerance=normalization_tolerance,
        )

    assert error_info.value.code == "categorical_probability_vector_tolerance_invalid"
    if math.isnan(normalization_tolerance):
        assert math.isnan(error_info.value.details["normalization_tolerance"])
    else:
        assert (
            error_info.value.details["normalization_tolerance"]
            == normalization_tolerance
        )
