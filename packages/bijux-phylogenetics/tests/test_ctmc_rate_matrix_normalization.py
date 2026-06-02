from __future__ import annotations

import numpy

from bijux_phylogenetics.phylo.likelihood import (
    compute_ctmc_expected_substitution_rate,
    normalize_ctmc_rate_matrix_by_expected_substitution_rate,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix


def test_normalize_ctmc_rate_matrix_by_expected_substitution_rate_is_scale_invariant() -> (
    None
):
    rate_matrix = numpy.array(
        [
            [-3.0, 2.0, 1.0],
            [4.0, -5.0, 1.0],
            [2.0, 3.0, -5.0],
        ],
        dtype=float,
    )
    stationary = numpy.array((11.0 / 21.0, 13.0 / 42.0, 1.0 / 6.0), dtype=float)

    normalized = normalize_ctmc_rate_matrix_by_expected_substitution_rate(
        rate_matrix,
        stationary,
        state_labels=("A", "B", "C"),
    )
    normalized_scaled = normalize_ctmc_rate_matrix_by_expected_substitution_rate(
        rate_matrix * 7.0,
        stationary,
        state_labels=("A", "B", "C"),
    )

    assert numpy.allclose(normalized, normalized_scaled, rtol=0.0, atol=1e-12)
    assert (
        compute_ctmc_expected_substitution_rate(
            normalized,
            stationary,
            state_labels=("A", "B", "C"),
        )
        == 1.0
    )


def test_transition_probability_matrix_matches_scaled_generator_after_branch_compensation() -> (
    None
):
    rate_matrix = numpy.array(
        [
            [-3.0, 2.0, 1.0],
            [4.0, -5.0, 1.0],
            [2.0, 3.0, -5.0],
        ],
        dtype=float,
    )
    stationary = numpy.array((11.0 / 21.0, 13.0 / 42.0, 1.0 / 6.0), dtype=float)
    branch_length = 0.375
    scale_factor = 11.0

    transition = transition_probability_matrix(rate_matrix, branch_length)
    scaled_transition = transition_probability_matrix(
        rate_matrix * scale_factor,
        branch_length / scale_factor,
    )
    normalized_transition = transition_probability_matrix(
        normalize_ctmc_rate_matrix_by_expected_substitution_rate(
            rate_matrix,
            stationary,
            state_labels=("A", "B", "C"),
        ),
        branch_length
        * compute_ctmc_expected_substitution_rate(
            rate_matrix,
            stationary,
            state_labels=("A", "B", "C"),
        ),
    )

    assert numpy.allclose(transition, scaled_transition, rtol=0.0, atol=1e-12)
    assert numpy.allclose(transition, normalized_transition, rtol=0.0, atol=1e-12)
