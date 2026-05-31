from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_f81_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_hky85_tree_likelihood_from_alignment,
    evaluate_k80_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_hky85_fixed_tree_likelihood_matches_independent_fixture() -> None:
    tree_path = fixture("trees", "hky85_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")

    report = evaluate_hky85_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        kappa=4.0,
        base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
    )
    expected_log_likelihood = _expected_hky85_two_tip_log_likelihood(
        base_frequencies=numpy.array([0.4, 0.1, 0.2, 0.3], dtype=float),
        kappa=4.0,
    )

    assert report.taxa == ["A", "B"]
    assert report.site_count == 4
    assert report.pattern_count == 4
    assert report.compression_used is True
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert report.base_frequency_source == "provided"
    assert math.isclose(report.base_frequency_a, 0.4, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.base_frequency_c, 0.1, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.base_frequency_g, 0.2, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.base_frequency_t, 0.3, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.kappa, 4.0, rel_tol=0.0, abs_tol=1e-12)
    assert report.parameter_count == 4
    assert math.isclose(
        report.log_likelihood,
        expected_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.aic,
        (-2.0 * expected_log_likelihood) + 8.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_hky85_reduces_to_k80_with_uniform_base_frequencies() -> None:
    tree = load_tree(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")
    )

    hky85_report = evaluate_hky85_tree_likelihood(
        tree,
        records,
        kappa=4.0,
        base_frequencies={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
    )
    k80_report = evaluate_k80_tree_likelihood(tree, records, kappa=4.0)

    assert math.isclose(
        hky85_report.log_likelihood,
        k80_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_hky85_reduces_to_f81_when_kappa_is_one() -> None:
    tree = load_tree(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")
    )
    base_frequencies = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}

    hky85_report = evaluate_hky85_tree_likelihood(
        tree,
        records,
        kappa=1.0,
        base_frequencies=base_frequencies,
    )
    f81_report = evaluate_f81_tree_likelihood(
        tree,
        records,
        base_frequencies=base_frequencies,
    )

    assert math.isclose(
        hky85_report.log_likelihood,
        f81_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _expected_hky85_two_tip_log_likelihood(
    *,
    base_frequencies: numpy.ndarray,
    kappa: float,
) -> float:
    rate_matrix = _independent_hky85_rate_matrix(
        base_frequencies=base_frequencies,
        kappa=kappa,
    )
    left_transition = transition_probability_matrix(rate_matrix, 0.1)
    right_transition = transition_probability_matrix(rate_matrix, 0.2)
    probability = 1.0
    state_index = {"A": 0, "C": 1, "G": 2, "T": 3}
    for left_state, right_state in (("A", "A"), ("C", "C"), ("A", "G"), ("A", "C")):
        pair_probability = 0.0
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        for root_index, root_probability in enumerate(base_frequencies):
            pair_probability += (
                float(root_probability)
                * float(left_transition[root_index, left_index])
                * float(right_transition[root_index, right_index])
            )
        probability *= pair_probability
    return math.log(probability)


def _independent_hky85_rate_matrix(
    *,
    base_frequencies: numpy.ndarray,
    kappa: float,
) -> numpy.ndarray:
    off_diagonal = numpy.zeros((4, 4), dtype=float)
    state_order = ("A", "C", "G", "T")
    transition_pairs = frozenset(
        {
            ("A", "G"),
            ("G", "A"),
            ("C", "T"),
            ("T", "C"),
        }
    )
    for left_index, left_state in enumerate(state_order):
        for right_index, right_state in enumerate(state_order):
            if left_index == right_index:
                continue
            multiplier = kappa if (left_state, right_state) in transition_pairs else 1.0
            off_diagonal[left_index, right_index] = (
                multiplier * base_frequencies[right_index]
            )
    rate_matrix = off_diagonal.copy()
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    expected_rate = -float(numpy.sum(base_frequencies * numpy.diag(rate_matrix)))
    return rate_matrix / expected_rate
