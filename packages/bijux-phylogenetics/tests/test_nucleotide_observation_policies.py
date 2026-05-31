from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_f81_tree_likelihood,
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
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
    assert math.isfinite(missing_report.log_likelihood)
    assert math.isfinite(ambiguity_report.log_likelihood)
    assert math.isfinite(fifth_state_report.log_likelihood)
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
