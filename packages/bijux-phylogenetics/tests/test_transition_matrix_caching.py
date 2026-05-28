from __future__ import annotations

import math

import numpy

from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
    transition_probability_matrix,
)


def test_transition_probability_matrix_returns_identity_at_zero_branch_length() -> None:
    matrix = transition_probability_matrix(
        numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float),
        0.0,
    )

    assert numpy.allclose(matrix, numpy.eye(2, dtype=float), rtol=0.0, atol=1e-12)


def test_transition_probability_matrix_matches_analytical_two_state_fixture() -> None:
    rate_matrix = numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float)
    branch_length = 0.25

    observed = transition_probability_matrix(rate_matrix, branch_length)

    assert numpy.allclose(
        observed,
        _analytical_two_state_transition_matrix(
            rate_a_to_b=2.0,
            rate_b_to_a=1.0,
            branch_length=branch_length,
        ),
        rtol=0.0,
        atol=1e-12,
    )


def test_transition_matrix_evaluator_reuses_cached_branch_lengths() -> None:
    rate_matrix = numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float)
    branch_lengths = [0.0, 0.25, 0.25, 0.5, 0.5, 0.0]

    cached_evaluator = build_transition_matrix_evaluator(
        rate_matrix,
        cache_matrices=True,
    )
    uncached_evaluator = build_transition_matrix_evaluator(
        rate_matrix,
        cache_matrices=False,
    )

    cached_matrices = [
        cached_evaluator.transition_probability_matrix(branch_length)
        for branch_length in branch_lengths
    ]
    uncached_matrices = [
        uncached_evaluator.transition_probability_matrix(branch_length)
        for branch_length in branch_lengths
    ]

    for cached_matrix, uncached_matrix in zip(
        cached_matrices,
        uncached_matrices,
        strict=True,
    ):
        assert numpy.allclose(cached_matrix, uncached_matrix, rtol=0.0, atol=1e-12)

    assert cached_evaluator.matrix_exponential_evaluation_count == 3
    assert cached_evaluator.cached_branch_length_count == 3
    assert uncached_evaluator.matrix_exponential_evaluation_count == len(branch_lengths)
    assert uncached_evaluator.cached_branch_length_count == 0


def test_transition_matrix_evaluator_keeps_model_specific_caches_separate() -> None:
    first_evaluator = build_transition_matrix_evaluator(
        numpy.array([[-2.0, 2.0], [1.0, -1.0]], dtype=float),
    )
    second_evaluator = build_transition_matrix_evaluator(
        numpy.array([[-1.0, 1.0], [0.25, -0.25]], dtype=float),
    )

    first_matrix = first_evaluator.transition_probability_matrix(0.25)
    second_matrix = second_evaluator.transition_probability_matrix(0.25)

    assert not numpy.allclose(first_matrix, second_matrix, rtol=0.0, atol=1e-12)
    assert first_evaluator.matrix_exponential_evaluation_count == 1
    assert second_evaluator.matrix_exponential_evaluation_count == 1


def _analytical_two_state_transition_matrix(
    *,
    rate_a_to_b: float,
    rate_b_to_a: float,
    branch_length: float,
) -> numpy.ndarray:
    total_rate = rate_a_to_b + rate_b_to_a
    decay = math.exp(-total_rate * branch_length)
    return numpy.array(
        [
            [
                (rate_b_to_a / total_rate) + (rate_a_to_b / total_rate) * decay,
                (rate_a_to_b / total_rate) - (rate_a_to_b / total_rate) * decay,
            ],
            [
                (rate_b_to_a / total_rate) - (rate_b_to_a / total_rate) * decay,
                (rate_a_to_b / total_rate) + (rate_b_to_a / total_rate) * decay,
            ],
        ],
        dtype=float,
    )
