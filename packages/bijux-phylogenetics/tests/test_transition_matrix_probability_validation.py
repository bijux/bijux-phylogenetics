from __future__ import annotations

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
    transition_probability_matrix,
)


def test_transition_probability_matrix_rejects_nonfinite_branch_lengths() -> None:
    rate_matrix = numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float)

    with pytest.raises(ValueError, match="branch length must be finite"):
        transition_probability_matrix(rate_matrix, float("inf"))

    with pytest.raises(ValueError, match="branch length must be finite"):
        transition_probability_matrix(rate_matrix, float("nan"))


def test_transition_matrix_evaluator_matches_direct_probability_across_scales() -> None:
    rate_matrix = numpy.array(
        [
            [-1.7, 0.8, 0.6, 0.3],
            [0.2, -1.1, 0.5, 0.4],
            [0.1, 0.3, -0.9, 0.5],
            [0.4, 0.2, 0.7, -1.3],
        ],
        dtype=float,
    )
    evaluator = build_transition_matrix_evaluator(rate_matrix)

    for branch_length in (1e-12, 0.25, 10.0, 1_000_000.0):
        direct = transition_probability_matrix(rate_matrix, branch_length)
        cached = evaluator.transition_probability_matrix(branch_length)

        assert numpy.allclose(cached, direct, rtol=0.0, atol=1e-12)
        assert numpy.all(cached >= 0.0)
        assert numpy.allclose(
            cached.sum(axis=1),
            numpy.ones(4, dtype=float),
            rtol=0.0,
            atol=1e-12,
        )
