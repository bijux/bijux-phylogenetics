from __future__ import annotations

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood import transition_probability_matrix
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_transition_probability_matrix_rejects_invalid_rate_matrix_before_exponentiation() -> (
    None
):
    with pytest.raises(PhylogeneticsError) as error_info:
        transition_probability_matrix(
            numpy.array([[-1.0, 0.9], [0.5, -0.5]], dtype=float),
            0.25,
        )

    assert error_info.value.code == "ctmc_rate_matrix_row_sums_nonzero"


def test_transition_probability_matrix_accepts_valid_ctmc_generator() -> None:
    matrix = transition_probability_matrix(
        numpy.array([[-1.0, 1.0], [0.25, -0.25]], dtype=float),
        0.25,
    )

    assert matrix.shape == (2, 2)
    assert numpy.allclose(matrix.sum(axis=1), numpy.ones(2, dtype=float))
