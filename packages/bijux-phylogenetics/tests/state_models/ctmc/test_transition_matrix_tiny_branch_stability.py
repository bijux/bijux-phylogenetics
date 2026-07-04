from __future__ import annotations

import numpy

from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix


def test_transition_probability_matrix_tiny_branch_tracks_first_order_generator() -> (
    None
):
    rate_matrix = numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float)
    branch_length = 1e-12

    observed = transition_probability_matrix(rate_matrix, branch_length)
    expected = numpy.eye(2, dtype=float) + (rate_matrix * branch_length)

    assert numpy.allclose(observed, expected, rtol=0.0, atol=1e-12)


def test_transition_probability_matrix_tiny_branch_stays_finite_on_asymmetric_four_state_generator() -> (
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
    branch_length = 1e-9

    observed = transition_probability_matrix(rate_matrix, branch_length)

    assert numpy.all(numpy.isfinite(observed))
    assert numpy.all(observed >= 0.0)
    assert numpy.allclose(
        observed.sum(axis=1),
        numpy.ones(4, dtype=float),
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.allclose(
        numpy.diag(observed),
        numpy.ones(4, dtype=float),
        rtol=0.0,
        atol=2e-9,
    )
