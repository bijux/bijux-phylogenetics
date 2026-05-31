from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_k80_tree_likelihood,
    evaluate_k80_tree_likelihood_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_k80_fixed_tree_likelihood_matches_two_tip_analytical_fixture() -> None:
    tree_path = fixture("trees", "k80_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", "k80_likelihood_alignment_2_taxa.fasta")

    report = evaluate_k80_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        kappa=4.0,
    )
    expected_probability = _expected_k80_two_tip_fixture_probability(kappa=4.0)

    assert report.taxa == ["A", "B"]
    assert report.site_count == 4
    assert report.pattern_count == 4
    assert report.compression_used is True
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert report.kappa == 4.0
    assert math.isclose(
        report.log_likelihood,
        math.log(expected_probability),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_k80_kappa_changes_likelihood_on_transition_biased_fixture() -> None:
    tree = load_tree(fixture("trees", "k80_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "k80_likelihood_alignment_2_taxa.fasta")
    )

    jc69_equivalent_report = evaluate_k80_tree_likelihood(tree, records, kappa=1.0)
    transition_weighted_report = evaluate_k80_tree_likelihood(tree, records, kappa=4.0)

    assert (
        transition_weighted_report.log_likelihood
        > jc69_equivalent_report.log_likelihood
    )


def _expected_k80_two_tip_fixture_probability(*, kappa: float) -> float:
    total_branch_length = 0.3
    beta = 1.0 / (kappa + 2.0)
    alpha = kappa * beta
    same_probability = 0.25 * (
        0.25
        + (0.25 * math.exp(-4.0 * beta * total_branch_length))
        + (0.5 * math.exp(-2.0 * (alpha + beta) * total_branch_length))
    )
    transition_probability = 0.25 * (
        0.25
        + (0.25 * math.exp(-4.0 * beta * total_branch_length))
        - (0.5 * math.exp(-2.0 * (alpha + beta) * total_branch_length))
    )
    transversion_probability = 0.25 * (
        0.25 - (0.25 * math.exp(-4.0 * beta * total_branch_length))
    )
    return (same_probability**2) * transition_probability * transversion_probability
