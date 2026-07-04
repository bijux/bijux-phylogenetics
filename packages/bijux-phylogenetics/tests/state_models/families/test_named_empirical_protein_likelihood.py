from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.likelihood import (
    EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
    evaluate_jc69_tree_likelihood_from_alignment,
    evaluate_named_empirical_protein_tree_likelihood_from_alignment,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_named_empirical_protein_models_reject_dna_surface_reuse() -> None:
    tree_path = fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
    )

    with pytest.raises(InvalidAlignmentError):
        evaluate_jc69_tree_likelihood_from_alignment(tree_path, alignment_path)

    report = evaluate_named_empirical_protein_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        model_name="LG",
    )

    assert report.matrix_label == "LG"
    assert report.root_prior_source == "stationary"
    assert report.state_count == 20
    assert (
        report.normalization_policy == EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY
    )
    assert math.isfinite(report.log_likelihood)


def test_named_empirical_protein_models_change_likelihood_on_shared_fixture() -> None:
    tree_path = fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
    )

    jtt_report = evaluate_named_empirical_protein_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        model_name="JTT",
    )
    wag_report = evaluate_named_empirical_protein_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        model_name="WAG",
    )
    lg_report = evaluate_named_empirical_protein_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        model_name="LG",
    )

    assert jtt_report.matrix_label == "JTT"
    assert wag_report.matrix_label == "WAG"
    assert lg_report.matrix_label == "LG"
    assert not math.isclose(
        jtt_report.log_likelihood,
        wag_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        wag_report.log_likelihood,
        lg_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        jtt_report.log_likelihood,
        lg_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
