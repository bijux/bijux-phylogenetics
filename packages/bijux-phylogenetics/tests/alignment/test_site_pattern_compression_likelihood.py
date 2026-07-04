from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    compress_alignment_site_patterns_from_records,
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
    sum_alignment_site_log_likelihoods,
    sum_compressed_site_pattern_log_likelihoods,
    transition_probability_matrix,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_compressed_jc69_site_patterns_match_uncompressed_likelihood() -> None:
    _assert_compressed_and_uncompressed_match(
        alignment_name="jc69_site_pattern_alignment.fasta",
        rate_matrix=_jc69_rate_matrix(),
        root_prior=numpy.full(4, 0.25, dtype=float),
    )


def test_compressed_hky_site_patterns_match_uncompressed_likelihood() -> None:
    _assert_compressed_and_uncompressed_match(
        alignment_name="hky_site_pattern_alignment.fasta",
        rate_matrix=_hky85_rate_matrix(
            base_frequencies=numpy.array([0.3, 0.2, 0.2, 0.3], dtype=float),
            kappa=4.0,
        ),
        root_prior=numpy.array([0.3, 0.2, 0.2, 0.3], dtype=float),
    )


def test_compressed_gtr_site_patterns_match_uncompressed_likelihood() -> None:
    _assert_compressed_and_uncompressed_match(
        alignment_name="gtr_site_pattern_alignment.fasta",
        rate_matrix=_gtr_rate_matrix(
            base_frequencies=numpy.array([0.28, 0.22, 0.24, 0.26], dtype=float),
            exchangeabilities={
                ("A", "C"): 0.8,
                ("A", "G"): 1.7,
                ("A", "T"): 0.5,
                ("C", "G"): 0.6,
                ("C", "T"): 1.2,
                ("G", "T"): 0.9,
            },
        ),
        root_prior=numpy.array([0.28, 0.22, 0.24, 0.26], dtype=float),
    )


def _assert_compressed_and_uncompressed_match(
    *,
    alignment_name: str,
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
) -> None:
    tree = load_tree(fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"))
    records = load_fasta_alignment(fixture("alignments", alignment_name))
    compressed = compress_alignment_site_patterns_from_records(records)
    taxon_order = [record.identifier for record in records]
    state_order = ["A", "C", "G", "T"]
    state_index = {state: index for index, state in enumerate(state_order)}
    transition_by_node_id = {
        child.node_id: transition_probability_matrix(
            rate_matrix,
            max(float(child.branch_length or 0.0), 0.0),
        )
        for _parent, child in tree.iter_edges()
    }

    def site_log_likelihood(states: tuple[str, ...]) -> float:
        states_by_taxon = dict(zip(taxon_order, states, strict=True))
        pruning_pass = postorder_conditional_likelihoods(
            tree,
            state_count=4,
            leaf_likelihood=lambda node: _one_hot(
                4,
                state_index[states_by_taxon[node.name or ""]],
            ),
            transition_matrix_for_child=lambda child: transition_by_node_id[
                child.node_id or ""
            ],
        )
        return log_likelihood_from_root_prior(
            tree,
            pruning_pass,
            root_prior=root_prior,
        )

    uncompressed = sum_alignment_site_log_likelihoods(
        records,
        site_log_likelihood=site_log_likelihood,
    )
    compressed_total = sum_compressed_site_pattern_log_likelihoods(
        compressed,
        site_log_likelihood=site_log_likelihood,
    )
    assert math.isclose(
        compressed_total,
        uncompressed,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _jc69_rate_matrix() -> numpy.ndarray:
    rate_matrix = numpy.full((4, 4), 1.0 / 3.0, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    return _finalize_rate_matrix(rate_matrix)


def _hky85_rate_matrix(
    *,
    base_frequencies: numpy.ndarray,
    kappa: float,
) -> numpy.ndarray:
    state_order = ["A", "C", "G", "T"]
    purines = {"A", "G"}
    pyrimidines = {"C", "T"}
    rate_matrix = numpy.zeros((4, 4), dtype=float)
    for left_index, left_state in enumerate(state_order):
        for right_index, right_state in enumerate(state_order):
            if left_index == right_index:
                continue
            is_transition = {left_state, right_state} <= purines or {
                left_state,
                right_state,
            } <= pyrimidines
            multiplier = kappa if is_transition else 1.0
            rate_matrix[left_index, right_index] = (
                multiplier * base_frequencies[right_index]
            )
    return _finalize_rate_matrix(rate_matrix)


def _gtr_rate_matrix(
    *,
    base_frequencies: numpy.ndarray,
    exchangeabilities: dict[tuple[str, str], float],
) -> numpy.ndarray:
    state_order = ["A", "C", "G", "T"]
    state_index = {state: index for index, state in enumerate(state_order)}
    rate_matrix = numpy.zeros((4, 4), dtype=float)
    for (left_state, right_state), exchangeability in exchangeabilities.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = (
            exchangeability * base_frequencies[right_index]
        )
        rate_matrix[right_index, left_index] = (
            exchangeability * base_frequencies[left_index]
        )
    return _finalize_rate_matrix(rate_matrix)


def _finalize_rate_matrix(off_diagonal_rates: numpy.ndarray) -> numpy.ndarray:
    rate_matrix = off_diagonal_rates.copy()
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    expected_rate = -float(
        numpy.sum(numpy.full(4, 0.25, dtype=float) * numpy.diag(rate_matrix))
    )
    return rate_matrix / expected_rate


def _one_hot(state_count: int, state_index: int) -> numpy.ndarray:
    vector = numpy.zeros(state_count, dtype=float)
    vector[state_index] = 1.0
    return vector
