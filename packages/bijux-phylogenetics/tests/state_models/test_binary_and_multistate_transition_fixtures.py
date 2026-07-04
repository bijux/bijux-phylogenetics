from __future__ import annotations

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood.ctmc import (
    solve_ctmc_stationary_distribution,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
    transition_probability_matrix,
)


@pytest.mark.parametrize(
    ("state_labels", "rate_matrix", "long_branch_length"),
    [
        (
            ("present", "absent"),
            numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float),
            25.0,
        ),
        (
            ("north", "south", "west", "island"),
            numpy.array(
                [
                    [-1.7, 0.8, 0.6, 0.3],
                    [0.2, -1.1, 0.5, 0.4],
                    [0.1, 0.3, -0.9, 0.5],
                    [0.4, 0.2, 0.7, -1.3],
                ],
                dtype=float,
            ),
            40.0,
        ),
    ],
)
def test_binary_and_multistate_transition_matrices_remain_valid_across_branch_scales(
    state_labels: tuple[str, ...],
    rate_matrix: numpy.ndarray,
    long_branch_length: float,
) -> None:
    tiny_branch_length = 1e-9
    evaluator = build_transition_matrix_evaluator(
        rate_matrix,
        cache_matrices=False,
    )
    stationary = numpy.array(
        solve_ctmc_stationary_distribution(
            rate_matrix,
            state_labels=state_labels,
        ).probabilities,
        dtype=float,
    )

    tiny_transition = evaluator.transition_probability_matrix(tiny_branch_length)
    long_transition = evaluator.transition_probability_matrix(long_branch_length)

    assert numpy.allclose(
        tiny_transition,
        numpy.eye(len(state_labels), dtype=float) + (rate_matrix * tiny_branch_length),
        rtol=0.0,
        atol=1e-9,
    )
    assert numpy.all(tiny_transition >= 0.0)
    assert numpy.all(long_transition >= 0.0)
    assert numpy.allclose(
        tiny_transition.sum(axis=1),
        numpy.ones(len(state_labels), dtype=float),
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.allclose(
        long_transition.sum(axis=1),
        numpy.ones(len(state_labels), dtype=float),
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.allclose(
        long_transition,
        numpy.repeat(stationary[numpy.newaxis, :], len(state_labels), axis=0),
        rtol=0.0,
        atol=1e-8,
    )


def test_binary_and_multistate_transition_evaluator_matches_direct_probability_surface() -> (
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
    evaluator = build_transition_matrix_evaluator(rate_matrix)

    for branch_length in (1e-12, 0.25, 7.5, 40.0):
        assert numpy.allclose(
            evaluator.transition_probability_matrix(branch_length),
            transition_probability_matrix(rate_matrix, branch_length),
            rtol=0.0,
            atol=1e-12,
        )
