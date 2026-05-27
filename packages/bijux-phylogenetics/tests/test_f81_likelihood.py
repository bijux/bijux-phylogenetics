from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_f81_tree_likelihood,
    evaluate_f81_tree_likelihood_from_alignment,
    evaluate_jc69_tree_likelihood,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_f81_fixed_tree_likelihood_matches_two_tip_analytical_fixture() -> None:
    tree_path = fixture("trees", "f81_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", "f81_likelihood_alignment_2_taxa.fasta")

    report = evaluate_f81_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
    )
    expected_log_likelihood = _expected_f81_two_tip_log_likelihood()

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
    assert report.parameter_count == 3
    assert math.isclose(
        report.log_likelihood,
        expected_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.aic,
        (-2.0 * expected_log_likelihood) + 6.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_f81_estimated_base_frequencies_prefer_biased_data_over_jc69_by_aic() -> None:
    tree = load_tree(fixture("trees", "f81_aic_bias_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "f81_aic_bias_alignment_2_taxa.fasta")
    )

    estimated_report = evaluate_f81_tree_likelihood(tree, records)
    uniform_report = evaluate_f81_tree_likelihood(
        tree,
        records,
        base_frequencies={"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
    )
    jc69_report = evaluate_jc69_tree_likelihood(tree, records)
    jc69_aic = -2.0 * jc69_report.log_likelihood

    assert estimated_report.base_frequency_source == "estimated"
    assert estimated_report.base_frequency_a > 0.25
    assert estimated_report.base_frequency_t < estimated_report.base_frequency_a
    assert estimated_report.base_frequency_c < 0.1
    assert estimated_report.base_frequency_g < 0.1
    assert estimated_report.aic < jc69_aic
    assert estimated_report.log_likelihood > jc69_report.log_likelihood
    assert math.isclose(
        uniform_report.log_likelihood,
        jc69_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert estimated_report.log_likelihood > uniform_report.log_likelihood


def _expected_f81_two_tip_log_likelihood() -> float:
    pi = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
    total_branch_length = 0.3
    variability = 1.0 - sum(value * value for value in pi.values())
    decay = math.exp(-total_branch_length / variability)
    probability = 1.0
    for left_state, right_state in (("A", "A"), ("C", "T"), ("T", "T"), ("G", "A")):
        if left_state == right_state:
            pair_probability = pi[left_state] * (
                pi[right_state] + ((1.0 - pi[right_state]) * decay)
            )
        else:
            pair_probability = pi[left_state] * (pi[right_state] * (1.0 - decay))
        probability *= pair_probability
    return math.log(probability)
