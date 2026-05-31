from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_nucleotide_marginal_ancestral_probabilities_from_alignment,
    jc69_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_STATE_INDEX,
    DNA_STATE_ORDER,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_marginal_ancestral_probabilities_match_jc69_two_tip_analytical_fixture() -> (
    None
):
    report = evaluate_nucleotide_marginal_ancestral_probabilities_from_alignment(
        fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta"),
        model_name="jc69",
    )

    assert report.model_name == "JC69"
    assert report.taxa == ["A", "B"]
    assert report.site_count == 4
    assert report.pattern_count == 4
    assert report.internal_node_count == 1
    assert report.compression_used is True
    assert report.expansion_policy == "expanded-internal-node-site-state-rows"
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert report.parameter_values == {}
    assert len(report.posterior_rows) == 16

    posterior_by_site_and_state = {
        (row.site_position, row.state): row.posterior_probability
        for row in report.posterior_rows
    }
    expected_by_site = {
        1: _expected_root_posterior(
            left_state="A",
            right_state="A",
            left_branch_length=0.1,
            right_branch_length=0.2,
        ),
        2: _expected_root_posterior(
            left_state="A",
            right_state="G",
            left_branch_length=0.1,
            right_branch_length=0.2,
        ),
        3: _expected_root_posterior(
            left_state="C",
            right_state="C",
            left_branch_length=0.1,
            right_branch_length=0.2,
        ),
        4: _expected_root_posterior(
            left_state="G",
            right_state="T",
            left_branch_length=0.1,
            right_branch_length=0.2,
        ),
    }
    for site_position, expected_vector in expected_by_site.items():
        for state_index, state in enumerate(DNA_STATE_ORDER):
            assert math.isclose(
                posterior_by_site_and_state[(site_position, state)],
                float(expected_vector[state_index]),
                rel_tol=0.0,
                abs_tol=1e-12,
            )


def test_marginal_ancestral_probabilities_expand_internal_node_site_state_rows() -> (
    None
):
    tree = load_tree(fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"))
    report = evaluate_nucleotide_marginal_ancestral_probabilities_from_alignment(
        fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
    )

    internal_node_ids = {
        node.node_id for node in tree.iter_internal_nodes(order="preorder")
    }
    assert internal_node_ids == {row.node_id for row in report.posterior_rows}
    assert report.internal_node_count == 3
    assert report.site_count == 10
    assert report.pattern_count == 6
    assert len(report.posterior_rows) == 120
    assert {row.site_position for row in report.posterior_rows} == set(range(1, 11))

    posterior_sums: dict[tuple[str, int], float] = {}
    for row in report.posterior_rows:
        key = (row.node_id, row.site_position)
        posterior_sums[key] = posterior_sums.get(key, 0.0) + row.posterior_probability
    assert len(posterior_sums) == report.internal_node_count * report.site_count
    for probability_sum in posterior_sums.values():
        assert math.isclose(probability_sum, 1.0, rel_tol=0.0, abs_tol=1e-12)


def test_marginal_ancestral_probabilities_normalize_for_all_selected_models() -> None:
    cases = [
        ("jc69", {}),
        ("k80", {"kappa": 2.0}),
        ("f81", {}),
        ("hky85", {"kappa": 2.5}),
        (
            "gtr",
            {
                "exchangeabilities": {
                    "AC": 1.0,
                    "AG": 2.0,
                    "AT": 1.5,
                    "CG": 0.8,
                    "CT": 1.7,
                    "GT": 1.2,
                }
            },
        ),
    ]

    for model_name, parameters in cases:
        report = evaluate_nucleotide_marginal_ancestral_probabilities_from_alignment(
            fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"),
            fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta"),
            model_name=model_name,
            **parameters,
        )
        assert len(report.posterior_rows) == 16
        assert {row.state for row in report.posterior_rows} == set(DNA_STATE_ORDER)
        posterior_sums: dict[tuple[str, int], float] = {}
        for row in report.posterior_rows:
            key = (row.node_id, row.site_position)
            posterior_sums[key] = (
                posterior_sums.get(key, 0.0) + row.posterior_probability
            )
        assert len(posterior_sums) == report.internal_node_count * report.site_count
        for probability_sum in posterior_sums.values():
            assert math.isclose(probability_sum, 1.0, rel_tol=0.0, abs_tol=1e-12)


def _expected_root_posterior(
    *,
    left_state: str,
    right_state: str,
    left_branch_length: float,
    right_branch_length: float,
) -> list[float]:
    weights = [0.0, 0.0, 0.0, 0.0]
    for state_index, state in enumerate(DNA_STATE_ORDER):
        left_probability = _jc69_transition_probability(
            start_state=state,
            end_state=left_state,
            branch_length=left_branch_length,
        )
        right_probability = _jc69_transition_probability(
            start_state=state,
            end_state=right_state,
            branch_length=right_branch_length,
        )
        weights[state_index] = 0.25 * left_probability * right_probability
    total = sum(weights)
    return [weight / total for weight in weights]


def _jc69_transition_probability(
    *,
    start_state: str,
    end_state: str,
    branch_length: float,
) -> float:
    matrix = jc69_transition_probability_matrix(branch_length)
    return float(matrix[DNA_STATE_INDEX[start_state], DNA_STATE_INDEX[end_state]])
