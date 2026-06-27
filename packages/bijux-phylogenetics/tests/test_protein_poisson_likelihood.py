from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_jc69_tree_likelihood_from_alignment,
    evaluate_protein_poisson_tree_likelihood,
    evaluate_protein_poisson_tree_likelihood_from_alignment,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_protein_poisson_fixed_tree_likelihood_matches_two_tip_analytical_fixture() -> (
    None
):
    tree_path = fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta"
    )

    report = evaluate_protein_poisson_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
    )
    expected_probability = _expected_two_tip_fixture_probability()

    assert report.taxa == ["A", "B"]
    assert report.site_count == 4
    assert report.pattern_count == 4
    assert report.compression_used is True
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert report.state_count == 20
    assert report.gap_policy == "treat-as-missing"
    assert report.missing_policy == "treat-as-missing"
    assert report.ambiguity_policy == "reject"
    assert report.root_prior_source == "stationary"
    assert math.isclose(
        report.log_likelihood,
        math.log(expected_probability),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_protein_poisson_gap_policy_can_reject_gap_states() -> None:
    tree_path = fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta"
    )

    with pytest.raises(InvalidAlignmentError, match="gap policy rejects '-'"):
        evaluate_protein_poisson_tree_likelihood_from_alignment(
            tree_path,
            alignment_path,
            gap_policy="reject",
        )


def test_protein_poisson_missing_policy_can_reject_missing_states() -> None:
    tree_path = fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta"
    )

    with pytest.raises(
        InvalidAlignmentError, match="missing-state policy rejects '\\?'"
    ):
        evaluate_protein_poisson_tree_likelihood_from_alignment(
            tree_path,
            alignment_path,
            missing_policy="reject",
        )


def test_protein_poisson_ambiguity_policy_can_reject_ambiguity_codes() -> None:
    tree_path = fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "protein_poisson_likelihood_alignment_ambiguity_2_taxa.fasta"
    )

    with pytest.raises(InvalidAlignmentError, match="ambiguity policy rejects 'B'"):
        evaluate_protein_poisson_tree_likelihood_from_alignment(
            tree_path,
            alignment_path,
            ambiguity_policy="reject",
        )


def test_protein_poisson_ambiguity_policies_change_likelihood() -> None:
    tree_path = fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "protein_poisson_likelihood_alignment_ambiguity_2_taxa.fasta"
    )

    missing_report = evaluate_protein_poisson_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        ambiguity_policy="treat-as-missing",
    )
    ambiguity_report = evaluate_protein_poisson_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        ambiguity_policy="ambiguity-vector",
    )

    assert missing_report.ambiguity_policy == "treat-as-missing"
    assert ambiguity_report.ambiguity_policy == "ambiguity-vector"
    assert math.isfinite(missing_report.log_likelihood)
    assert math.isfinite(ambiguity_report.log_likelihood)
    assert not math.isclose(
        missing_report.log_likelihood,
        ambiguity_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_protein_poisson_likelihood_does_not_reuse_dna_alphabet_surface() -> None:
    tree_path = fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta"
    )

    with pytest.raises(InvalidAlignmentError):
        evaluate_jc69_tree_likelihood_from_alignment(tree_path, alignment_path)

    report = evaluate_protein_poisson_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
    )

    assert math.isfinite(report.log_likelihood)


def test_protein_poisson_root_prior_policies_change_asymmetric_likelihoods() -> None:
    tree_path = fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta"
    )

    equal_report = evaluate_protein_poisson_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        root_prior_policy="equal",
    )
    empirical_report = evaluate_protein_poisson_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        root_prior_policy="empirical",
    )
    stationary_report = evaluate_protein_poisson_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        root_prior_policy="stationary",
    )
    fixed_report = evaluate_protein_poisson_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        root_prior_policy="fixed-state",
        fixed_root_state="A",
    )
    provided_report = evaluate_protein_poisson_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        root_prior_policy="provided",
        root_prior=_biased_protein_root_prior(),
    )

    assert equal_report.root_prior_source == "equal"
    assert empirical_report.root_prior_source == "empirical"
    assert stationary_report.root_prior_source == "stationary"
    assert fixed_report.root_prior_source == "fixed-state"
    assert provided_report.root_prior_source == "provided"
    assert math.isclose(
        equal_report.log_likelihood,
        stationary_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        empirical_report.log_likelihood,
        equal_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        fixed_report.log_likelihood,
        equal_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        provided_report.log_likelihood,
        equal_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _biased_protein_root_prior() -> numpy.ndarray:
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


def _expected_two_tip_fixture_probability() -> float:
    total_branch_length = 0.3
    state_count = 20.0
    decay = math.exp((-state_count * total_branch_length) / (state_count - 1.0))
    same_site_probability = (1.0 / state_count) * (
        (1.0 / state_count) + (((state_count - 1.0) / state_count) * decay)
    )
    different_site_probability = (1.0 / state_count) * (
        (1.0 / state_count) - ((1.0 / state_count) * decay)
    )
    single_observed_with_missing_probability = 1.0 / state_count
    return (
        same_site_probability
        * different_site_probability
        * single_observed_with_missing_probability
        * single_observed_with_missing_probability
    )
