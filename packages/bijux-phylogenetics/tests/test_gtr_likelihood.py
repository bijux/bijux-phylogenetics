from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_gtr_tree_likelihood,
    evaluate_gtr_tree_likelihood_from_alignment,
    evaluate_hky85_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_gtr_fixed_tree_likelihood_matches_independent_fixture() -> None:
    tree_path = fixture("trees", "gtr_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", "gtr_likelihood_alignment_2_taxa.fasta")
    exchangeabilities = {
        "AC": 1.0,
        "AG": 4.5,
        "AT": 0.8,
        "CG": 1.6,
        "CT": 2.4,
        "GT": 3.1,
    }
    base_frequencies = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}

    report = evaluate_gtr_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        exchangeabilities=exchangeabilities,
        base_frequencies=base_frequencies,
    )
    expected_log_likelihood = _expected_gtr_two_tip_log_likelihood(
        exchangeabilities=numpy.array([1.0, 4.5, 0.8, 1.6, 2.4, 3.1], dtype=float),
        base_frequencies=numpy.array([0.4, 0.1, 0.2, 0.3], dtype=float),
    )

    assert report.taxa == ["A", "B"]
    assert report.site_count == 4
    assert report.pattern_count == 4
    assert report.compression_used is True
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert report.base_frequency_source == "provided"
    assert report.exchangeability_anchor == "AC=1"
    assert math.isclose(report.exchangeability_ac, 1.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.exchangeability_ag, 4.5, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.exchangeability_at, 0.8, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.exchangeability_cg, 1.6, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.exchangeability_ct, 2.4, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.exchangeability_gt, 3.1, rel_tol=0.0, abs_tol=1e-12)
    assert report.parameter_count == 8
    assert math.isclose(
        report.log_likelihood,
        expected_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.aic,
        (-2.0 * expected_log_likelihood) + 16.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_gtr_reduces_to_hky85_with_hky_shaped_exchangeabilities() -> None:
    tree = load_tree(fixture("trees", "gtr_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "gtr_likelihood_alignment_2_taxa.fasta")
    )
    base_frequencies = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
    hky_shaped_exchangeabilities = {
        "AC": 1.0,
        "AG": 4.0,
        "AT": 1.0,
        "CG": 1.0,
        "CT": 4.0,
        "GT": 1.0,
    }

    gtr_report = evaluate_gtr_tree_likelihood(
        tree,
        records,
        exchangeabilities=hky_shaped_exchangeabilities,
        base_frequencies=base_frequencies,
    )
    hky_report = evaluate_hky85_tree_likelihood(
        tree,
        records,
        kappa=4.0,
        base_frequencies=base_frequencies,
    )

    assert math.isclose(
        gtr_report.log_likelihood,
        hky_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _expected_gtr_two_tip_log_likelihood(
    *,
    exchangeabilities: numpy.ndarray,
    base_frequencies: numpy.ndarray,
) -> float:
    rate_matrix = _independent_gtr_rate_matrix(
        exchangeabilities=exchangeabilities,
        base_frequencies=base_frequencies,
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


def _independent_gtr_rate_matrix(
    *,
    exchangeabilities: numpy.ndarray,
    base_frequencies: numpy.ndarray,
) -> numpy.ndarray:
    rate_matrix = numpy.zeros((4, 4), dtype=float)
    state_order = ("A", "C", "G", "T")
    pair_order = (
        ("A", "C"),
        ("A", "G"),
        ("A", "T"),
        ("C", "G"),
        ("C", "T"),
        ("G", "T"),
    )
    for pair_index, (left_state, right_state) in enumerate(pair_order):
        left_index = state_order.index(left_state)
        right_index = state_order.index(right_state)
        exchangeability = float(exchangeabilities[pair_index])
        rate_matrix[left_index, right_index] = exchangeability * float(
            base_frequencies[right_index]
        )
        rate_matrix[right_index, left_index] = exchangeability * float(
            base_frequencies[left_index]
        )
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    expected_rate = -float(numpy.sum(base_frequencies * numpy.diag(rate_matrix)))
    return rate_matrix / expected_rate
