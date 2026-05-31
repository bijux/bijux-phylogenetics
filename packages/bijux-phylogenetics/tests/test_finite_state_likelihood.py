from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.ancestral.discrete.likelihood import (
    transition_probability_matrix as discrete_transition_probability_matrix,
)
from bijux_phylogenetics.ancestral.discrete.likelihood import (
    tree_log_likelihood as discrete_tree_log_likelihood,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


def test_two_tip_finite_state_pruning_matches_binary_analytical_likelihood() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    transition_by_taxon = {
        "A": numpy.array([[0.9, 0.1], [0.2, 0.8]], dtype=float),
        "B": numpy.array([[0.9, 0.1], [0.2, 0.8]], dtype=float),
    }
    observed_state_index = {"A": 0, "B": 1}

    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=2,
        leaf_likelihood=lambda node: _one_hot(
            2,
            observed_state_index[node.name or ""],
        ),
        transition_matrix_for_child=lambda child: transition_by_taxon[child.name or ""],
    )
    log_likelihood = log_likelihood_from_root_prior(
        tree,
        pruning_pass,
        root_prior=numpy.array([0.5, 0.5], dtype=float),
    )

    assert numpy.allclose(
        pruning_pass.conditional_for_node(tree.root),
        numpy.array([0.36, 0.64], dtype=float),
        rtol=0.0,
        atol=1e-12,
    )
    assert math.isclose(
        log_likelihood,
        math.log(0.125),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_three_tip_finite_state_pruning_matches_dna_like_analytical_likelihood() -> (
    None
):
    tree = load_tree(fixture("felsenstein_three_tip_tree.nwk"))
    transition_by_descendants = {
        ("A",): numpy.array(
            [
                [0.7, 0.1, 0.1, 0.1],
                [0.2, 0.5, 0.2, 0.1],
                [0.1, 0.2, 0.6, 0.1],
                [0.25, 0.25, 0.25, 0.25],
            ],
            dtype=float,
        ),
        ("B", "C"): numpy.array(
            [
                [0.6, 0.2, 0.1, 0.1],
                [0.1, 0.7, 0.1, 0.1],
                [0.15, 0.15, 0.6, 0.1],
                [0.25, 0.25, 0.25, 0.25],
            ],
            dtype=float,
        ),
        ("B",): numpy.array(
            [
                [0.65, 0.15, 0.1, 0.1],
                [0.2, 0.5, 0.2, 0.1],
                [0.1, 0.25, 0.55, 0.1],
                [0.25, 0.25, 0.25, 0.25],
            ],
            dtype=float,
        ),
        ("C",): numpy.array(
            [
                [0.55, 0.15, 0.2, 0.1],
                [0.1, 0.6, 0.2, 0.1],
                [0.15, 0.15, 0.6, 0.1],
                [0.25, 0.25, 0.25, 0.25],
            ],
            dtype=float,
        ),
    }
    observed_state_index = {"A": 0, "B": 1, "C": 2}

    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=4,
        leaf_likelihood=lambda node: _one_hot(
            4,
            observed_state_index[node.name or ""],
        ),
        transition_matrix_for_child=lambda child: transition_by_descendants[
            tuple(child.descendant_taxa)
        ],
    )
    log_likelihood = log_likelihood_from_root_prior(
        tree,
        pruning_pass,
        root_prior=numpy.full(4, 0.25, dtype=float),
    )

    assert math.isclose(
        math.exp(log_likelihood),
        0.0233265625,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_finite_state_pruning_handles_twenty_state_protein_like_kernel() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    state_count = 20
    transition = numpy.zeros((state_count, state_count), dtype=float)
    for index in range(state_count):
        transition[index, index] = 0.55
        transition[index, (index + 1) % state_count] = 0.2
        transition[index, (index + 5) % state_count] = 0.15
        transition[index, (index + 9) % state_count] = 0.1
    observed_state_index = {"A": 3, "B": 8}

    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=state_count,
        leaf_likelihood=lambda node: _one_hot(
            state_count,
            observed_state_index[node.name or ""],
        ),
        transition_matrix_for_child=lambda _child: transition,
    )
    log_likelihood = log_likelihood_from_root_prior(
        tree,
        pruning_pass,
        root_prior=numpy.full(state_count, 1.0 / state_count, dtype=float),
    )

    assert pruning_pass.state_count == 20
    assert pruning_pass.conditional_for_node(tree.root).shape == (20,)
    assert math.isclose(
        math.exp(log_likelihood),
        0.004125,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_discrete_tree_log_likelihood_matches_shared_pruning_kernel() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    states_by_taxon = {"A": "0", "B": "1"}
    rate_matrix = numpy.array([[-1.0, 1.0], [2.0, -2.0]], dtype=float)
    state_order = ["0", "1"]
    root_prior = numpy.array([0.5, 0.5], dtype=float)

    shared_pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=2,
        leaf_likelihood=lambda node: _one_hot(
            2,
            0 if states_by_taxon[node.name or ""] == "0" else 1,
        ),
        transition_matrix_for_child=lambda child: (
            discrete_transition_probability_matrix(
                rate_matrix,
                max(float(child.branch_length or 0.0), 0.0),
            )
        ),
    )
    shared_log_likelihood = log_likelihood_from_root_prior(
        tree,
        shared_pruning_pass,
        root_prior=root_prior,
    )
    discrete_log_likelihood = discrete_tree_log_likelihood(
        tree,
        states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )

    assert math.isclose(
        discrete_log_likelihood,
        shared_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _one_hot(state_count: int, state_index: int) -> numpy.ndarray:
    vector = numpy.zeros(state_count, dtype=float)
    vector[state_index] = 1.0
    return vector
