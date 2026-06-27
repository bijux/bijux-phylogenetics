from __future__ import annotations

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood.codon_states import (
    build_equal_rate_codon_ctmc_rate_matrix,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    solve_ctmc_stationary_distribution,
)
from bijux_phylogenetics.phylo.likelihood.poisson import (
    protein_poisson_rate_matrix,
    protein_poisson_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
    transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.protein import PROTEIN_STATE_ORDER


@pytest.mark.parametrize(
    ("state_labels", "rate_matrix", "long_branch_length", "long_branch_tolerance"),
    [
        (
            PROTEIN_STATE_ORDER,
            protein_poisson_rate_matrix(),
            20.0,
            1e-8,
        ),
        (
            build_equal_rate_codon_ctmc_rate_matrix()[0].state_order,
            build_equal_rate_codon_ctmc_rate_matrix()[1],
            40.0,
            2e-8,
        ),
    ],
)
def test_protein_and_codon_transition_matrices_remain_valid_across_branch_scales(
    state_labels: tuple[str, ...],
    rate_matrix: numpy.ndarray,
    long_branch_length: float,
    long_branch_tolerance: float,
) -> None:
    tiny_branch_length = 1e-9
    evaluator = build_transition_matrix_evaluator(rate_matrix)
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
        atol=long_branch_tolerance,
    )


def test_protein_and_codon_transition_evaluator_matches_direct_probability_surface() -> (
    None
):
    state_space, codon_rate_matrix, _root_prior, _frequency_source = (
        build_equal_rate_codon_ctmc_rate_matrix()
    )
    evaluator = build_transition_matrix_evaluator(codon_rate_matrix)

    for branch_length in (1e-12, 0.25, 5.0, 15.0):
        direct = transition_probability_matrix(codon_rate_matrix, branch_length)
        cached = evaluator.transition_probability_matrix(branch_length)

        assert direct.shape == (
            len(state_space.state_order),
            len(state_space.state_order),
        )
        assert numpy.allclose(cached, direct, rtol=0.0, atol=1e-12)


def test_protein_poisson_transition_matches_normalized_ctmc_engine() -> None:
    rate_matrix = protein_poisson_rate_matrix()

    for branch_length in (0.0, 1e-12, 0.25, 5.0, 20.0):
        native = protein_poisson_transition_probability_matrix(branch_length)
        generic = transition_probability_matrix(rate_matrix, branch_length)

        assert native.shape == (20, 20)
        assert numpy.all(native >= 0.0)
        assert numpy.allclose(
            native.sum(axis=1),
            numpy.ones(20, dtype=float),
            rtol=0.0,
            atol=1e-12,
        )
        assert numpy.allclose(native, generic, rtol=0.0, atol=1e-12)
