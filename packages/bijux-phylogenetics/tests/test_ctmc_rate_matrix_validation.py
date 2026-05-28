from __future__ import annotations

import math

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood.ctmc import (
    validate_ctmc_rate_matrix,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_validate_ctmc_rate_matrix_preserves_valid_asymmetric_matrix() -> None:
    validated = validate_ctmc_rate_matrix(
        [
            [-0.3, 0.2, 0.1],
            [0.4, -0.5, 0.1],
            [0.2, 0.0, -0.2],
        ],
        state_labels=("A", "B", "C"),
    )

    assert validated.state_labels == ("A", "B", "C")
    assert validated.state_count == 3
    assert validated.index_for_state("B") == 1
    assert numpy.allclose(
        validated.rate_matrix,
        numpy.array(
            [
                [-0.3, 0.2, 0.1],
                [0.4, -0.5, 0.1],
                [0.2, 0.0, -0.2],
            ],
            dtype=float,
        ),
    )


def test_validate_ctmc_rate_matrix_assigns_default_state_labels() -> None:
    validated = validate_ctmc_rate_matrix([[-1.0, 1.0], [0.5, -0.5]])

    assert validated.state_labels == ("state_0", "state_1")


def test_validate_ctmc_rate_matrix_rejects_non_square_input() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        validate_ctmc_rate_matrix([[0.0, 1.0]])

    assert error_info.value.code == "ctmc_rate_matrix_not_square"
    assert error_info.value.details["shape"] == [1, 2]


@pytest.mark.parametrize("invalid_value", (math.nan, math.inf, -math.inf))
def test_validate_ctmc_rate_matrix_rejects_non_finite_values(
    invalid_value: float,
) -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        validate_ctmc_rate_matrix([[-1.0, 1.0], [invalid_value, -0.5]])

    assert error_info.value.code == "ctmc_rate_matrix_value_not_finite"
    assert error_info.value.details["row_index"] == 1
    assert error_info.value.details["column_index"] == 0


def test_validate_ctmc_rate_matrix_rejects_negative_off_diagonal_rate() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        validate_ctmc_rate_matrix([[-1.0, 1.0], [-0.1, 0.1]])

    assert error_info.value.code == "ctmc_rate_matrix_off_diagonal_negative"
    assert error_info.value.details["source_state"] == "state_1"
    assert error_info.value.details["target_state"] == "state_0"
    assert error_info.value.details["value"] == -0.1


def test_validate_ctmc_rate_matrix_rejects_nonzero_row_sums() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        validate_ctmc_rate_matrix(
            [[-0.8, 0.7, 0.1], [0.2, -0.5, 0.2], [0.0, 0.0, 0.0]],
            state_labels=("A", "B", "C"),
        )

    assert error_info.value.code == "ctmc_rate_matrix_row_sums_nonzero"
    assert error_info.value.details["offending_rows"] == [
        {"row_index": 1, "state_label": "B", "row_sum": -0.09999999999999998}
    ]


def test_validate_ctmc_rate_matrix_rejects_duplicate_state_labels() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        validate_ctmc_rate_matrix(
            [[-1.0, 1.0], [0.5, -0.5]],
            state_labels=("A", "A"),
        )

    assert error_info.value.code == "ctmc_rate_matrix_state_labels_duplicate"
    assert error_info.value.details["duplicate_state_labels"] == ["A"]


def test_validate_ctmc_rate_matrix_rejects_state_label_count_mismatch() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        validate_ctmc_rate_matrix(
            [[-1.0, 1.0], [0.5, -0.5]],
            state_labels=("A", "B", "C"),
        )

    assert error_info.value.code == "ctmc_rate_matrix_state_label_count_mismatch"
    assert error_info.value.details["state_count"] == 2
    assert error_info.value.details["state_label_count"] == 3
