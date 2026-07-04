from __future__ import annotations

import numpy

from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix


def test_transition_probability_matrix_long_two_state_branch_converges_to_stationary_rows() -> (
    None
):
    rate_matrix = numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float)

    observed = transition_probability_matrix(rate_matrix, 1_000_000.0)
    expected_stationary_row = numpy.array([1.0 / 3.0, 2.0 / 3.0], dtype=float)

    assert numpy.allclose(
        observed[0],
        expected_stationary_row,
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.allclose(
        observed[1],
        expected_stationary_row,
        rtol=0.0,
        atol=1e-12,
    )


def test_transition_probability_matrix_long_asymmetric_branch_preserves_valid_rows() -> (
    None
):
    rate_matrix = numpy.array(
        [
            [-1.7, 0.8, 0.6, 0.3],
            [0.2, -1.1, 0.5, 0.4],
            [0.1, 0.3, -0.9, 0.5],
            [0.4, 0.2, 0.7, -1.3],
        ],
        dtype=float,
    )

    observed = transition_probability_matrix(rate_matrix, 1_000_000.0)

    assert numpy.all(numpy.isfinite(observed))
    assert numpy.all(observed >= 0.0)
    assert numpy.allclose(
        observed.sum(axis=1),
        numpy.ones(4, dtype=float),
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.allclose(
        observed[0],
        observed[1],
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.allclose(
        observed[1],
        observed[2],
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.allclose(
        observed[2],
        observed[3],
        rtol=0.0,
        atol=1e-12,
    )
