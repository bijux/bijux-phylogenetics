from __future__ import annotations

import math

import numpy

from bijux_phylogenetics.phylo.likelihood import (
    build_equal_rate_codon_ctmc_rate_matrix,
    compute_ctmc_expected_substitution_rate,
)
from bijux_phylogenetics.phylo.likelihood.dna import (
    validate_dna_base_frequencies,
)
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    augment_dna_rate_matrix_with_gap_state,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix


def test_augment_dna_rate_matrix_with_gap_state_is_scale_invariant() -> None:
    stationary = validate_dna_base_frequencies(
        {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
        model_name="GTR",
    )
    nucleotide_rate_matrix = numpy.array(
        [
            [-1.1, 0.1, 0.4, 0.6],
            [0.4, -0.9, 0.2, 0.3],
            [0.8, 0.1, -1.2, 0.3],
            [0.8, 0.1, 0.2, -1.1],
        ],
        dtype=float,
    )

    augmented = augment_dna_rate_matrix_with_gap_state(
        nucleotide_rate_matrix,
        nucleotide_frequencies=stationary,
        gap_state_frequency=0.15,
        gap_exchangeability=0.75,
        model_name="GTR",
    )
    augmented_scaled = augment_dna_rate_matrix_with_gap_state(
        nucleotide_rate_matrix * 8.0,
        nucleotide_frequencies=stationary,
        gap_state_frequency=0.15,
        gap_exchangeability=0.75 * 8.0,
        model_name="GTR",
    )

    assert numpy.allclose(augmented, augmented_scaled, rtol=0.0, atol=1e-12)


def test_equal_rate_codon_ctmc_rate_matrix_preserves_scale_compensated_transitions() -> (
    None
):
    state_space, rate_matrix, frequencies, _ = build_equal_rate_codon_ctmc_rate_matrix()
    scale_factor = 13.0
    branch_length = 0.7

    assert math.isclose(
        compute_ctmc_expected_substitution_rate(
            rate_matrix,
            frequencies,
            state_labels=state_space.state_order,
        ),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert numpy.allclose(
        transition_probability_matrix(rate_matrix, branch_length),
        transition_probability_matrix(
            rate_matrix * scale_factor,
            branch_length / scale_factor,
        ),
        rtol=0.0,
        atol=1e-12,
    )
