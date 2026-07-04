from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood.dna import (
    evaluate_fixed_topology_dna_site_log_likelihood,
    normalize_dna_rate_matrix,
    normalize_unambiguous_dna_records,
    validate_dna_base_frequencies,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix

FIXTURES = Path(__file__).parent / "fixtures"
_F81_BASE_FREQUENCIES = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_normalize_dna_rate_matrix_is_scale_invariant() -> None:
    stationary = validate_dna_base_frequencies(_F81_BASE_FREQUENCIES, model_name="F81")
    off_diagonal_rates = numpy.array(
        [
            [0.0, 0.1, 0.2, 0.3],
            [0.4, 0.0, 0.2, 0.1],
            [0.3, 0.1, 0.0, 0.2],
            [0.2, 0.3, 0.1, 0.0],
        ],
        dtype=float,
    )

    normalized = normalize_dna_rate_matrix(
        off_diagonal_rates,
        stationary_frequencies=stationary,
        model_name="F81",
    )
    normalized_scaled = normalize_dna_rate_matrix(
        off_diagonal_rates * 9.0,
        stationary_frequencies=stationary,
        model_name="F81",
    )

    assert numpy.allclose(normalized, normalized_scaled, rtol=0.0, atol=1e-12)


def test_f81_site_log_likelihood_is_invariant_under_rate_scaling_and_branch_compensation() -> (
    None
):
    tree = load_tree(fixture("trees", "f81_likelihood_tree_2_taxa.nwk"))
    records = normalize_unambiguous_dna_records(
        load_fasta_alignment(
            fixture("alignments", "f81_likelihood_alignment_2_taxa.fasta")
        ),
        model_name="F81 normalization test",
    )
    taxon_order = [record.identifier for record in records]
    stationary = validate_dna_base_frequencies(_F81_BASE_FREQUENCIES, model_name="F81")
    off_diagonal_rates = numpy.array(
        [
            [0.0, stationary[1], stationary[2], stationary[3]],
            [stationary[0], 0.0, stationary[2], stationary[3]],
            [stationary[0], stationary[1], 0.0, stationary[3]],
            [stationary[0], stationary[1], stationary[2], 0.0],
        ],
        dtype=float,
    )
    normalized_rate_matrix = normalize_dna_rate_matrix(
        off_diagonal_rates,
        stationary_frequencies=stationary,
        model_name="F81",
    )
    scaled_rate_matrix = normalized_rate_matrix * 5.0

    def total_log_likelihood(
        rate_matrix: numpy.ndarray, *, branch_scale: float
    ) -> float:
        site_count = len(records[0].sequence)
        total = 0.0
        for site_index in range(site_count):
            states = tuple(record.sequence[site_index] for record in records)
            total += evaluate_fixed_topology_dna_site_log_likelihood(
                tree,
                states,
                taxon_order=taxon_order,
                model_name="F81 normalization test",
                observation_policy="reject",
                root_prior=stationary,
                transition_matrix_for_child=lambda child: transition_probability_matrix(
                    rate_matrix,
                    float(child.branch_length or 0.0) / branch_scale,
                ),
            )
        return total

    original = total_log_likelihood(normalized_rate_matrix, branch_scale=1.0)
    scaled = total_log_likelihood(scaled_rate_matrix, branch_scale=5.0)

    assert math.isclose(original, scaled, rel_tol=0.0, abs_tol=1e-12)
