from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    compute_marginal_state_posteriors,
    jc69_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_STATE_INDEX,
    DNA_STATE_ORDER,
    UNIFORM_DNA_ROOT_PRIOR,
    one_hot_dna_leaf_vector,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_finite_state_posteriors_match_jc69_two_tip_same_state_fixture() -> None:
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    states_by_taxon = {"A": "A", "B": "A"}

    posterior_pass = compute_marginal_state_posteriors(
        tree,
        state_count=4,
        leaf_likelihood=lambda node: one_hot_dna_leaf_vector(
            states_by_taxon,
            model_name="JC69",
            node_name=node.name,
        ),
        transition_matrix_for_child=lambda child: jc69_transition_probability_matrix(
            max(float(child.branch_length or 0.0), 0.0)
        ),
        root_prior=UNIFORM_DNA_ROOT_PRIOR,
    )

    expected = _expected_root_posterior(
        left_state="A",
        right_state="A",
        left_branch_length=0.1,
        right_branch_length=0.2,
    )
    observed = posterior_pass.posterior_for_node(tree.root)

    assert math.isclose(float(observed.sum()), 1.0, rel_tol=0.0, abs_tol=1e-12)
    assert numpy.allclose(observed, expected, rtol=0.0, atol=1e-12)


def test_finite_state_posteriors_match_jc69_two_tip_different_state_fixture() -> None:
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    states_by_taxon = {"A": "A", "B": "G"}

    posterior_pass = compute_marginal_state_posteriors(
        tree,
        state_count=4,
        leaf_likelihood=lambda node: one_hot_dna_leaf_vector(
            states_by_taxon,
            model_name="JC69",
            node_name=node.name,
        ),
        transition_matrix_for_child=lambda child: jc69_transition_probability_matrix(
            max(float(child.branch_length or 0.0), 0.0)
        ),
        root_prior=UNIFORM_DNA_ROOT_PRIOR,
    )

    expected = _expected_root_posterior(
        left_state="A",
        right_state="G",
        left_branch_length=0.1,
        right_branch_length=0.2,
    )
    observed = posterior_pass.posterior_for_node(tree.root)

    assert math.isclose(float(observed.sum()), 1.0, rel_tol=0.0, abs_tol=1e-12)
    assert numpy.allclose(observed, expected, rtol=0.0, atol=1e-12)


def _expected_root_posterior(
    *,
    left_state: str,
    right_state: str,
    left_branch_length: float,
    right_branch_length: float,
) -> numpy.ndarray:
    weights = numpy.zeros(4, dtype=float)
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
    return weights / float(weights.sum())


def _jc69_transition_probability(
    *,
    start_state: str,
    end_state: str,
    branch_length: float,
) -> float:
    matrix = jc69_transition_probability_matrix(branch_length)
    return float(matrix[DNA_STATE_INDEX[start_state], DNA_STATE_INDEX[end_state]])
