from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    compress_alignment_site_patterns_from_records,
    sum_alignment_site_log_likelihoods,
    sum_compressed_site_pattern_log_likelihoods,
    transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.protein import (
    evaluate_fixed_topology_protein_site_log_likelihood,
    normalize_unambiguous_protein_records,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_compressed_protein_patterns_match_uncompressed_likelihood() -> None:
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = normalize_unambiguous_protein_records(
        load_fasta_alignment(
            fixture(
                "alignments",
                "empirical_protein_site_pattern_alignment_2_taxa.fasta",
            )
        ),
        model_name="empirical protein matrix",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    rate_matrix = _compact_polar_rate_matrix()
    root_prior = _biased_root_prior()
    transition_by_node_id = {
        child.node_id: transition_probability_matrix(
            rate_matrix,
            max(float(child.branch_length or 0.0), 0.0),
        )
        for _parent, child in tree.iter_edges()
    }

    def site_log_likelihood(states: tuple[str, ...]) -> float:
        return evaluate_fixed_topology_protein_site_log_likelihood(
            tree,
            states,
            taxon_order=compressed_patterns.taxon_order,
            model_name="empirical protein matrix",
            root_prior=root_prior,
            transition_matrix_for_child=lambda child: transition_by_node_id[
                child.node_id or ""
            ],
        )

    uncompressed_total = sum_alignment_site_log_likelihoods(
        records,
        site_log_likelihood=site_log_likelihood,
    )
    compressed_total = sum_compressed_site_pattern_log_likelihoods(
        compressed_patterns,
        site_log_likelihood=site_log_likelihood,
    )

    assert compressed_patterns.pattern_count == 3
    assert [pattern.weight for pattern in compressed_patterns.patterns] == [2, 3, 1]
    assert math.isclose(
        compressed_total,
        uncompressed_total,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _compact_polar_rate_matrix() -> numpy.ndarray:
    state_order = _protein_state_order()
    state_index = {state: index for index, state in enumerate(state_order)}
    rate_matrix = numpy.full((len(state_order), len(state_order)), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for (left_state, right_state), rate in {
        ("A", "C"): 0.45,
        ("C", "D"): 0.35,
        ("D", "E"): 0.55,
        ("A", "E"): 0.20,
    }.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix


def _biased_root_prior() -> numpy.ndarray:
    prior = numpy.full(20, 0.02, dtype=float)
    state_index = {state: index for index, state in enumerate(_protein_state_order())}
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())


def _protein_state_order() -> tuple[str, ...]:
    return (
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    )
