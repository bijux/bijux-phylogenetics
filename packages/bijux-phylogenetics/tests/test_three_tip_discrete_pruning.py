from __future__ import annotations

from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    build_discrete_mk_finite_state_surface,
    compile_finite_state_likelihood_context,
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from tests.support.three_tip_pruning_reference import (
    build_three_tip_pruning_reference,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


@pytest.mark.parametrize(
    ("model_name", "state_order", "states_by_taxon", "rate_matrix"),
    [
        (
            "unordered discrete Mk",
            ["north", "south", "west"],
            {"A": "north", "B": "south", "C": "west"},
            rate_matrix_from_log_parameters(
                numpy.log(numpy.array([0.8], dtype=float)),
                state_order=["north", "south", "west"],
                model="equal-rates",
                state_ordering="unordered",
                allowed_transition_pairs={
                    (0, 1),
                    (0, 2),
                    (1, 0),
                    (1, 2),
                    (2, 0),
                    (2, 1),
                },
            ),
        ),
        (
            "ordered discrete Mk",
            ["north", "south", "island"],
            {"A": "north", "B": "south", "C": "island"},
            rate_matrix_from_log_parameters(
                numpy.log(numpy.array([0.6], dtype=float)),
                state_order=["north", "south", "island"],
                model="equal-rates",
                state_ordering="ordered",
                allowed_transition_pairs={
                    (0, 1),
                    (1, 0),
                    (1, 2),
                    (2, 1),
                },
            ),
        ),
        (
            "irreversible discrete Mk",
            ["present", "absent"],
            {"A": "present", "B": "present", "C": "absent"},
            rate_matrix_from_log_parameters(
                numpy.log(numpy.array([0.4], dtype=float)),
                state_order=["present", "absent"],
                model="irreversible",
                state_ordering="unordered",
                allowed_transition_pairs={(0, 1)},
            ),
        ),
    ],
)
def test_discrete_surfaces_match_three_tip_reference(
    model_name: str,
    state_order: list[str],
    states_by_taxon: dict[str, str],
    rate_matrix: numpy.ndarray,
) -> None:
    root_prior = numpy.full(len(state_order), 1.0 / len(state_order), dtype=float)
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        model_name=model_name,
    )

    assert numpy.allclose(surface.root_prior, root_prior, rtol=0.0, atol=1e-12)
    _assert_surface_matches_three_tip_reference(
        surface=surface,
        states_by_taxon=states_by_taxon,
    )


def _assert_surface_matches_three_tip_reference(
    *,
    surface,
    states_by_taxon: dict[str, str],
) -> None:
    tree = load_tree(fixture("felsenstein_three_tip_tree.nwk"))
    context = compile_finite_state_likelihood_context(tree, surface)
    leaf_vector_by_taxon = {
        taxon: surface.leaf_likelihood_vector_for_state(taxon, state)
        for taxon, state in states_by_taxon.items()
    }
    reference = build_three_tip_pruning_reference(
        tree=tree,
        root_prior=surface.root_prior,
        leaf_vector_by_taxon=leaf_vector_by_taxon,
        transition_matrix_by_node_id=context.transition_matrix_by_node_id,
    )
    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=surface.state_count,
        leaf_likelihood=lambda node: leaf_vector_by_taxon[node.name or ""],
        transition_matrix_for_child=lambda child: context.transition_matrix_by_node_id[
            child.node_id
        ],
    )
    log_likelihood = log_likelihood_from_root_prior(
        tree,
        pruning_pass,
        root_prior=surface.root_prior,
    )

    assert numpy.allclose(
        pruning_pass.conditional_for_node(reference.internal_child),
        reference.internal_normalized,
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.allclose(
        pruning_pass.conditional_for_node(tree.root),
        reference.root_normalized,
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.isclose(
        pruning_pass.subtree_log_scaling_for_node(reference.internal_child),
        reference.internal_log_scale,
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.isclose(
        pruning_pass.subtree_log_scaling_for_node(tree.root),
        reference.root_log_scale,
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.isclose(log_likelihood, reference.log_likelihood, rtol=0.0, atol=1e-12)
