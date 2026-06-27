from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.codon import (
    evaluate_codon_ctmc_tree_likelihood,
    evaluate_codon_ctmc_tree_likelihood_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.codon_observation_policies import (
    resolve_codon_observation_leaf_vector,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import (
    build_equal_rate_codon_ctmc_rate_matrix,
    resolve_codon_state_space,
    validate_codon_frequency_vector,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from tests.support.two_tip_likelihood_reference import (
    two_tip_log_likelihood_for_labels,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_codon_ctmc_likelihood_matches_two_tip_high_precision_fixture() -> None:
    report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta"),
    )
    expected_log_likelihood = _expected_codon_two_tip_log_likelihood()

    assert report.taxa == ["A", "B"]
    assert report.site_count == 2
    assert report.pattern_count == 2
    assert report.compression_used is False
    assert report.state_count == 61
    assert report.genetic_code_id == 1
    assert report.genetic_code_name == "Standard"
    assert report.codon_frequency_source == "uniform"
    assert report.root_prior_source == "stationary"
    assert report.observation_policy == "reject"
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert math.isclose(
        report.log_likelihood,
        expected_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_standard_codon_state_space_excludes_stop_codons_and_keeps_61_states() -> None:
    state_space = resolve_codon_state_space()

    assert len(state_space.state_order) == 61
    assert state_space.stop_codons == frozenset({"TAA", "TAG", "TGA"})
    assert "TAA" not in state_space.state_order
    assert "TAG" not in state_space.state_order
    assert "TGA" not in state_space.state_order
    assert state_space.state_order[0] == "AAA"
    assert state_space.state_order[-1] == "TTT"


def test_codon_ctmc_likelihood_rejects_stop_codon_states() -> None:
    try:
        evaluate_codon_ctmc_tree_likelihood_from_alignment(
            fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
            fixture("alignments", "codon_likelihood_alignment_stop_2_taxa.fasta"),
        )
    except InvalidAlignmentError as error:
        assert "excludes stop codon states" in str(error)
        assert "TAA" in str(error)
    else:
        raise AssertionError("stop-codon codon site should be rejected")


def test_codon_ctmc_likelihood_changes_with_codon_frequencies() -> None:
    uniform_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta"),
    )
    state_space = resolve_codon_state_space()
    biased_frequencies = dict.fromkeys(state_space.state_order, 1.0)
    biased_frequencies["AAA"] = 20.0
    biased_frequencies["AAG"] = 15.0
    biased_frequencies["GCT"] = 0.5
    biased_frequencies["GCC"] = 0.5

    biased_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta"),
        codon_frequencies=biased_frequencies,
    )

    assert biased_report.codon_frequency_source == "provided"
    assert not math.isclose(
        biased_report.log_likelihood,
        uniform_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_codon_frequency_vector_rejects_stop_codon_keys() -> None:
    state_space = resolve_codon_state_space()
    invalid_frequencies = dict.fromkeys(state_space.state_order, 1.0)
    invalid_frequencies["TAA"] = 1.0

    with pytest.raises(
        InvalidAlignmentError,
        match="unexpected codons: TAA",
    ):
        validate_codon_frequency_vector(
            invalid_frequencies,
            state_space=state_space,
            owner_name="codon CTMC likelihood",
        )


def test_codon_ctmc_root_prior_policies_change_asymmetric_likelihoods() -> None:
    tree_path = fixture("trees", "codon_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    state_space = resolve_codon_state_space()
    biased_frequencies = dict.fromkeys(state_space.state_order, 1.0)
    biased_frequencies["AAA"] = 20.0
    biased_frequencies["AAG"] = 15.0
    biased_frequencies["GCT"] = 0.5
    biased_frequencies["GCC"] = 0.5
    provided_root_prior = dict.fromkeys(state_space.state_order, 1.0)
    provided_root_prior["AAA"] = 25.0
    provided_root_prior["AAG"] = 4.0
    provided_root_prior["GCT"] = 0.2
    provided_root_prior["GCC"] = 0.2

    equal_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        codon_frequencies=biased_frequencies,
        root_prior_policy="equal",
    )
    empirical_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        codon_frequencies=biased_frequencies,
        root_prior_policy="empirical",
    )
    stationary_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        codon_frequencies=biased_frequencies,
        root_prior_policy="stationary",
    )
    fixed_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        codon_frequencies=biased_frequencies,
        root_prior_policy="fixed-state",
        fixed_root_state="AAA",
    )
    provided_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        codon_frequencies=biased_frequencies,
        root_prior_policy="provided",
        root_prior=provided_root_prior,
    )

    assert equal_report.root_prior_source == "equal"
    assert empirical_report.root_prior_source == "empirical"
    assert stationary_report.root_prior_source == "stationary"
    assert fixed_report.root_prior_source == "fixed-state"
    assert provided_report.root_prior_source == "provided"
    assert not math.isclose(
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


def test_codon_ctmc_likelihood_rejects_ambiguous_codons_under_reject_policy() -> None:
    with pytest.raises(
        InvalidAlignmentError,
        match="observation policy 'reject' requires resolved sense codons only",
    ):
        evaluate_codon_ctmc_tree_likelihood_from_alignment(
            fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
            fixture("alignments", "codon_likelihood_alignment_ambiguity_2_taxa.fasta"),
            observation_policy="reject",
        )


def test_codon_ctmc_observation_policies_change_likelihood_on_ambiguity_fixture() -> (
    None
):
    missing_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_ambiguity_2_taxa.fasta"),
        observation_policy="treat-as-missing",
    )
    ambiguity_report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_ambiguity_2_taxa.fasta"),
        observation_policy="ambiguity-vector",
    )

    assert missing_report.observation_policy == "treat-as-missing"
    assert ambiguity_report.observation_policy == "ambiguity-vector"
    assert math.isfinite(missing_report.log_likelihood)
    assert math.isfinite(ambiguity_report.log_likelihood)
    assert not math.isclose(
        missing_report.log_likelihood,
        ambiguity_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_codon_ambiguity_vectors_preserve_partial_missing_and_gap_constraints() -> None:
    state_space = resolve_codon_state_space()

    partial_missing_vector = resolve_codon_observation_leaf_vector(
        "AT?",
        model_name="codon CTMC",
        node_name="A",
        observation_policy="ambiguity-vector",
        state_order=state_space.state_order,
        state_index=state_space.state_index,
    )
    partial_gap_vector = resolve_codon_observation_leaf_vector(
        "A-G",
        model_name="codon CTMC",
        node_name="A",
        observation_policy="ambiguity-vector",
        state_order=state_space.state_order,
        state_index=state_space.state_index,
    )
    treat_as_missing_vector = resolve_codon_observation_leaf_vector(
        "AT?",
        model_name="codon CTMC",
        node_name="A",
        observation_policy="treat-as-missing",
        state_order=state_space.state_order,
        state_index=state_space.state_index,
    )

    expected_partial_missing = {"ATA", "ATC", "ATG", "ATT"}
    expected_partial_gap = {"AAG", "ACG", "AGG", "ATG"}

    assert {
        codon
        for codon, index in state_space.state_index.items()
        if partial_missing_vector[index] == 1.0
    } == expected_partial_missing
    assert {
        codon
        for codon, index in state_space.state_index.items()
        if partial_gap_vector[index] == 1.0
    } == expected_partial_gap
    assert partial_missing_vector.sum() == 4.0
    assert partial_gap_vector.sum() == 4.0
    assert treat_as_missing_vector.sum() == float(len(state_space.state_order))


def test_codon_partial_missing_and_gap_change_likelihood_under_ambiguity_vectors() -> (
    None
):
    records = [
        AlignmentRecord(identifier="A", sequence="AT?A-G"),
        AlignmentRecord(identifier="B", sequence="ATGAAG"),
    ]

    missing_report = evaluate_codon_ctmc_tree_likelihood(
        load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk")),
        records,
        observation_policy="treat-as-missing",
    )
    ambiguity_report = evaluate_codon_ctmc_tree_likelihood(
        load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk")),
        records,
        observation_policy="ambiguity-vector",
    )

    assert math.isfinite(missing_report.log_likelihood)
    assert math.isfinite(ambiguity_report.log_likelihood)
    assert not math.isclose(
        missing_report.log_likelihood,
        ambiguity_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_codon_ctmc_observation_policy_treats_missing_and_gap_codons_explicitly() -> (
    None
):
    records = [
        AlignmentRecord(identifier="A", sequence="ATG???GGG"),
        AlignmentRecord(identifier="B", sequence="ATG---GGG"),
    ]

    with pytest.raises(
        InvalidAlignmentError,
        match="does not allow missing or gap codon state",
    ):
        evaluate_codon_ctmc_tree_likelihood(
            load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk")),
            records,
            observation_policy="reject",
        )

    report = evaluate_codon_ctmc_tree_likelihood(
        load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk")),
        records,
        observation_policy="treat-as-missing",
    )

    assert report.observation_policy == "treat-as-missing"
    assert math.isfinite(report.log_likelihood)


def _expected_codon_two_tip_log_likelihood() -> float:
    state_space, rate_matrix, root_prior, _frequency_source = (
        build_equal_rate_codon_ctmc_rate_matrix()
    )
    left_transition = transition_probability_matrix(rate_matrix, 0.1)
    right_transition = transition_probability_matrix(rate_matrix, 0.2)
    return two_tip_log_likelihood_for_labels(
        observed_pairs=(("AAA", "AAG"), ("GCT", "GCC")),
        state_index=state_space.state_index,
        root_prior=root_prior,
        left_transition=left_transition,
        right_transition=right_transition,
    )
