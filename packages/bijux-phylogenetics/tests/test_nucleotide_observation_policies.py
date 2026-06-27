from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_f81_tree_likelihood,
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    estimate_empirical_dna_base_frequencies_from_records,
    estimate_empirical_gap_state_frequency,
    resolve_dna_observation_leaf_vector,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_jc69_observation_policies_change_likelihood_on_ambiguity_fixture() -> None:
    tree = load_tree(
        fixture("trees", "jc69_joint_ancestral_difference_tree_3_taxa.nwk")
    )
    records = load_fasta_alignment(
        fixture("alignments", "example_alignment_ambiguity.fasta")
    )

    missing_report = evaluate_jc69_tree_likelihood(
        tree,
        records,
        observation_policy="treat-as-missing",
    )
    ambiguity_report = evaluate_jc69_tree_likelihood(
        tree,
        records,
        observation_policy="ambiguity-vector",
    )
    fifth_state_report = evaluate_jc69_tree_likelihood(
        tree,
        records,
        observation_policy="fifth-state",
    )

    assert missing_report.state_count == 4
    assert missing_report.observation_policy == "treat-as-missing"
    assert ambiguity_report.state_count == 4
    assert ambiguity_report.observation_policy == "ambiguity-vector"
    assert fifth_state_report.state_count == 5
    assert fifth_state_report.observation_policy == "fifth-state"
    assert missing_report.parameter_count == 0
    assert ambiguity_report.parameter_count == 0
    assert fifth_state_report.parameter_count == 0
    assert math.isfinite(missing_report.log_likelihood)
    assert math.isfinite(ambiguity_report.log_likelihood)
    assert math.isfinite(fifth_state_report.log_likelihood)
    assert math.isclose(
        missing_report.aic,
        -2.0 * missing_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        ambiguity_report.aic,
        -2.0 * ambiguity_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        fifth_state_report.aic,
        -2.0 * fifth_state_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        missing_report.log_likelihood,
        ambiguity_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        ambiguity_report.log_likelihood,
        fifth_state_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_k80_observation_policies_preserve_one_parameter_likelihood_contract() -> (
    None
):
    tree = load_tree(
        fixture("trees", "jc69_joint_ancestral_difference_tree_3_taxa.nwk")
    )
    records = load_fasta_alignment(
        fixture("alignments", "example_alignment_ambiguity.fasta")
    )

    missing_report = evaluate_k80_tree_likelihood(
        tree,
        records,
        kappa=3.0,
        observation_policy="treat-as-missing",
    )
    ambiguity_report = evaluate_k80_tree_likelihood(
        tree,
        records,
        kappa=3.0,
        observation_policy="ambiguity-vector",
    )
    fifth_state_report = evaluate_k80_tree_likelihood(
        tree,
        records,
        kappa=3.0,
        observation_policy="fifth-state",
    )

    assert missing_report.state_count == 4
    assert missing_report.observation_policy == "treat-as-missing"
    assert ambiguity_report.state_count == 4
    assert ambiguity_report.observation_policy == "ambiguity-vector"
    assert fifth_state_report.state_count == 5
    assert fifth_state_report.observation_policy == "fifth-state"
    assert missing_report.parameter_count == 1
    assert ambiguity_report.parameter_count == 1
    assert fifth_state_report.parameter_count == 1
    assert math.isfinite(missing_report.log_likelihood)
    assert math.isfinite(ambiguity_report.log_likelihood)
    assert math.isfinite(fifth_state_report.log_likelihood)
    assert math.isclose(
        missing_report.aic,
        (-2.0 * missing_report.log_likelihood) + 2.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        ambiguity_report.aic,
        (-2.0 * ambiguity_report.log_likelihood) + 2.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        fifth_state_report.aic,
        (-2.0 * fifth_state_report.log_likelihood) + 2.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        missing_report.log_likelihood,
        ambiguity_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        ambiguity_report.log_likelihood,
        fifth_state_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_f81_observation_policies_preserve_frequency_parameter_contract() -> None:
    tree = load_tree(
        fixture("trees", "jc69_joint_ancestral_difference_tree_3_taxa.nwk")
    )
    records = load_fasta_alignment(
        fixture("alignments", "example_alignment_ambiguity.fasta")
    )

    missing_report = evaluate_f81_tree_likelihood(
        tree,
        records,
        observation_policy="treat-as-missing",
    )
    ambiguity_report = evaluate_f81_tree_likelihood(
        tree,
        records,
        observation_policy="ambiguity-vector",
    )
    fifth_state_report = evaluate_f81_tree_likelihood(
        tree,
        records,
        observation_policy="fifth-state",
    )

    assert missing_report.state_count == 4
    assert missing_report.observation_policy == "treat-as-missing"
    assert ambiguity_report.state_count == 4
    assert ambiguity_report.observation_policy == "ambiguity-vector"
    assert fifth_state_report.state_count == 5
    assert fifth_state_report.observation_policy == "fifth-state"
    assert missing_report.parameter_count == 3
    assert ambiguity_report.parameter_count == 3
    assert fifth_state_report.parameter_count == 3
    assert missing_report.base_frequency_source == "estimated"
    assert ambiguity_report.base_frequency_source == "estimated"
    assert fifth_state_report.base_frequency_source == "estimated"
    assert math.isfinite(missing_report.log_likelihood)
    assert math.isfinite(ambiguity_report.log_likelihood)
    assert math.isfinite(fifth_state_report.log_likelihood)
    assert math.isclose(
        missing_report.aic,
        (-2.0 * missing_report.log_likelihood) + 6.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        ambiguity_report.aic,
        (-2.0 * ambiguity_report.log_likelihood) + 6.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        fifth_state_report.aic,
        (-2.0 * fifth_state_report.log_likelihood) + 6.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        missing_report.log_likelihood,
        ambiguity_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        ambiguity_report.log_likelihood,
        fifth_state_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_hky85_observation_policies_preserve_full_parameter_contract() -> None:
    tree = load_tree(
        fixture("trees", "jc69_joint_ancestral_difference_tree_3_taxa.nwk")
    )
    records = load_fasta_alignment(
        fixture("alignments", "example_alignment_ambiguity.fasta")
    )

    missing_report = evaluate_hky85_tree_likelihood(
        tree,
        records,
        kappa=3.0,
        observation_policy="treat-as-missing",
    )
    ambiguity_report = evaluate_hky85_tree_likelihood(
        tree,
        records,
        kappa=3.0,
        observation_policy="ambiguity-vector",
    )
    fifth_state_report = evaluate_hky85_tree_likelihood(
        tree,
        records,
        kappa=3.0,
        observation_policy="fifth-state",
    )

    assert missing_report.state_count == 4
    assert missing_report.observation_policy == "treat-as-missing"
    assert ambiguity_report.state_count == 4
    assert ambiguity_report.observation_policy == "ambiguity-vector"
    assert fifth_state_report.state_count == 5
    assert fifth_state_report.observation_policy == "fifth-state"
    assert missing_report.parameter_count == 4
    assert ambiguity_report.parameter_count == 4
    assert fifth_state_report.parameter_count == 4
    assert missing_report.base_frequency_source == "estimated"
    assert ambiguity_report.base_frequency_source == "estimated"
    assert fifth_state_report.base_frequency_source == "estimated"
    assert math.isfinite(missing_report.log_likelihood)
    assert math.isfinite(ambiguity_report.log_likelihood)
    assert math.isfinite(fifth_state_report.log_likelihood)
    assert math.isclose(
        missing_report.aic,
        (-2.0 * missing_report.log_likelihood) + 8.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        ambiguity_report.aic,
        (-2.0 * ambiguity_report.log_likelihood) + 8.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        fifth_state_report.aic,
        (-2.0 * fifth_state_report.log_likelihood) + 8.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        missing_report.log_likelihood,
        ambiguity_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        ambiguity_report.log_likelihood,
        fifth_state_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_gtr_observation_policies_preserve_full_parameter_contract() -> None:
    tree = load_tree(
        fixture("trees", "jc69_joint_ancestral_difference_tree_3_taxa.nwk")
    )
    records = load_fasta_alignment(
        fixture("alignments", "example_alignment_ambiguity.fasta")
    )
    exchangeabilities = {
        "AC": 1.0,
        "AG": 2.0,
        "AT": 0.5,
        "CG": 1.5,
        "CT": 1.75,
        "GT": 1.25,
    }

    missing_report = evaluate_gtr_tree_likelihood(
        tree,
        records,
        exchangeabilities=exchangeabilities,
        observation_policy="treat-as-missing",
    )
    ambiguity_report = evaluate_gtr_tree_likelihood(
        tree,
        records,
        exchangeabilities=exchangeabilities,
        observation_policy="ambiguity-vector",
    )
    fifth_state_report = evaluate_gtr_tree_likelihood(
        tree,
        records,
        exchangeabilities=exchangeabilities,
        observation_policy="fifth-state",
    )

    assert missing_report.state_count == 4
    assert missing_report.observation_policy == "treat-as-missing"
    assert ambiguity_report.state_count == 4
    assert ambiguity_report.observation_policy == "ambiguity-vector"
    assert fifth_state_report.state_count == 5
    assert fifth_state_report.observation_policy == "fifth-state"
    assert missing_report.parameter_count == 8
    assert ambiguity_report.parameter_count == 8
    assert fifth_state_report.parameter_count == 8
    assert missing_report.base_frequency_source == "estimated"
    assert ambiguity_report.base_frequency_source == "estimated"
    assert fifth_state_report.base_frequency_source == "estimated"
    assert math.isfinite(missing_report.log_likelihood)
    assert math.isfinite(ambiguity_report.log_likelihood)
    assert math.isfinite(fifth_state_report.log_likelihood)
    assert math.isclose(
        missing_report.aic,
        (-2.0 * missing_report.log_likelihood) + 16.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        ambiguity_report.aic,
        (-2.0 * ambiguity_report.log_likelihood) + 16.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        fifth_state_report.aic,
        (-2.0 * fifth_state_report.log_likelihood) + 16.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        missing_report.log_likelihood,
        ambiguity_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        ambiguity_report.log_likelihood,
        fifth_state_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_dna_observation_vectors_remain_explicit_across_policies() -> None:
    ambiguity_vector = resolve_dna_observation_leaf_vector(
        {"A": "R"},
        model_name="jc69",
        node_name="A",
        observation_policy="ambiguity-vector",
    )
    ambiguity_as_missing_vector = resolve_dna_observation_leaf_vector(
        {"A": "R"},
        model_name="jc69",
        node_name="A",
        observation_policy="treat-as-missing",
    )
    fifth_state_gap_vector = resolve_dna_observation_leaf_vector(
        {"A": "-"},
        model_name="jc69",
        node_name="A",
        observation_policy="fifth-state",
    )
    missing_vector = resolve_dna_observation_leaf_vector(
        {"A": "?"},
        model_name="jc69",
        node_name="A",
        observation_policy="fifth-state",
    )

    numpy.testing.assert_array_equal(ambiguity_vector, numpy.array([1, 0, 1, 0]))
    numpy.testing.assert_array_equal(
        ambiguity_as_missing_vector,
        numpy.ones(4, dtype=float),
    )
    numpy.testing.assert_array_equal(
        fifth_state_gap_vector,
        numpy.array([0, 0, 0, 0, 1]),
    )
    numpy.testing.assert_array_equal(missing_vector, numpy.ones(5, dtype=float))


def test_dna_empirical_base_frequencies_weight_ambiguity_vectors_explicitly() -> None:
    records = [
        AlignmentRecord(identifier="A", sequence="R"),
        AlignmentRecord(identifier="B", sequence="C"),
        AlignmentRecord(identifier="C", sequence="T"),
    ]

    ambiguity_weighted = estimate_empirical_dna_base_frequencies_from_records(
        records,
        model_name="jc69",
        observation_policy="ambiguity-vector",
    )
    missing_weighted = estimate_empirical_dna_base_frequencies_from_records(
        records,
        model_name="jc69",
        observation_policy="treat-as-missing",
    )

    numpy.testing.assert_allclose(
        ambiguity_weighted,
        numpy.array([1.0 / 6.0, 1.0 / 3.0, 1.0 / 6.0, 1.0 / 3.0]),
    )
    numpy.testing.assert_allclose(
        missing_weighted,
        numpy.array([0.0, 0.5, 0.0, 0.5]),
    )


def test_dna_gap_frequency_estimation_counts_only_informative_fifth_state_sites() -> (
    None
):
    records = [
        AlignmentRecord(identifier="A", sequence="A-?"),
        AlignmentRecord(identifier="B", sequence="CR?"),
    ]

    gap_frequency = estimate_empirical_gap_state_frequency(
        records,
        model_name="jc69",
    )

    assert math.isclose(gap_frequency, 0.25, rel_tol=0.0, abs_tol=1e-12)


@pytest.mark.parametrize(
    ("evaluator", "kwargs"),
    [
        (evaluate_jc69_tree_likelihood, {}),
        (evaluate_k80_tree_likelihood, {"kappa": 3.0}),
        (evaluate_f81_tree_likelihood, {}),
        (evaluate_hky85_tree_likelihood, {"kappa": 3.0}),
        (
            evaluate_gtr_tree_likelihood,
            {
                "exchangeabilities": {
                    "AC": 1.0,
                    "AG": 2.0,
                    "AT": 0.5,
                    "CG": 1.5,
                    "CT": 1.75,
                    "GT": 1.25,
                }
            },
        ),
    ],
)
def test_nucleotide_models_support_fifth_state_observation_policy(
    evaluator,
    kwargs: dict[str, object],
) -> None:
    tree = load_tree(
        fixture("trees", "jc69_joint_ancestral_difference_tree_3_taxa.nwk")
    )
    records = load_fasta_alignment(
        fixture("alignments", "example_alignment_ambiguity.fasta")
    )

    report = evaluator(
        tree,
        records,
        observation_policy="fifth-state",
        **kwargs,
    )

    assert report.state_count == 5
    assert report.observation_policy == "fifth-state"
    assert math.isfinite(report.log_likelihood)
    if evaluator is evaluate_jc69_tree_likelihood:
        assert report.parameter_count == 0
        assert math.isclose(
            report.aic,
            -2.0 * report.log_likelihood,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    if evaluator is evaluate_k80_tree_likelihood:
        assert report.parameter_count == 1
        assert math.isclose(
            report.aic,
            (-2.0 * report.log_likelihood) + 2.0,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    if evaluator is evaluate_f81_tree_likelihood:
        assert report.parameter_count == 3
        assert math.isclose(
            report.aic,
            (-2.0 * report.log_likelihood) + 6.0,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    if evaluator is evaluate_hky85_tree_likelihood:
        assert report.parameter_count == 4
        assert math.isclose(
            report.aic,
            (-2.0 * report.log_likelihood) + 8.0,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    if evaluator is evaluate_gtr_tree_likelihood:
        assert report.parameter_count == 8
        assert math.isclose(
            report.aic,
            (-2.0 * report.log_likelihood) + 16.0,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
