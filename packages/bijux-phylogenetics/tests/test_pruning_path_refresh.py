from __future__ import annotations

import math

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood.pruning import (
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
    refresh_postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def test_refresh_postorder_conditional_likelihoods_updates_only_affected_path() -> (
    None
):
    tree = _balanced_four_tip_tree()
    states_by_taxon = {"A": "0", "B": "1", "C": "0", "D": "1"}
    transition_by_node_id = {
        node.node_id: _transition_matrix(float(node.branch_length or 0.0))
        for node in tree.iter_nodes(order="preorder")
        if node is not tree.root
    }

    original_pass = postorder_conditional_likelihoods(
        tree,
        state_count=2,
        leaf_likelihood=lambda node: _leaf_vector(states_by_taxon[node.name]),
        transition_matrix_for_child=lambda child: transition_by_node_id[child.node_id],
    )
    original_root_log_likelihood = log_likelihood_from_root_prior(
        tree,
        original_pass,
        root_prior=numpy.array([0.5, 0.5], dtype=float),
    )

    left_internal = tree.root.children[0]
    right_internal = tree.root.children[1]
    changed_leaf = left_internal.children[0]
    transition_by_node_id[changed_leaf.node_id] = _transition_matrix(0.35)

    unaffected_conditional = original_pass.conditional_by_node_id[right_internal.node_id]
    unaffected_scale = original_pass.subtree_log_scaling_by_node_id[right_internal.node_id]
    previous_left_conditional = original_pass.conditional_by_node_id[left_internal.node_id]
    previous_root_conditional = original_pass.conditional_by_node_id[tree.root.node_id]

    refresh_postorder_conditional_likelihoods(
        tree,
        original_pass,
        affected_node_ids=(left_internal.node_id, tree.root.node_id),
        transition_matrix_for_child=lambda child: transition_by_node_id[child.node_id],
    )

    refreshed_root_log_likelihood = log_likelihood_from_root_prior(
        tree,
        original_pass,
        root_prior=numpy.array([0.5, 0.5], dtype=float),
    )
    recomputed_pass = postorder_conditional_likelihoods(
        tree,
        state_count=2,
        leaf_likelihood=lambda node: _leaf_vector(states_by_taxon[node.name]),
        transition_matrix_for_child=lambda child: transition_by_node_id[child.node_id],
    )
    recomputed_root_log_likelihood = log_likelihood_from_root_prior(
        tree,
        recomputed_pass,
        root_prior=numpy.array([0.5, 0.5], dtype=float),
    )

    assert not math.isclose(
        original_root_log_likelihood,
        refreshed_root_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        refreshed_root_log_likelihood,
        recomputed_root_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert original_pass.conditional_by_node_id[right_internal.node_id] is unaffected_conditional
    assert (
        original_pass.subtree_log_scaling_by_node_id[right_internal.node_id]
        == unaffected_scale
    )
    assert (
        original_pass.conditional_by_node_id[left_internal.node_id]
        is not previous_left_conditional
    )
    assert (
        original_pass.conditional_by_node_id[tree.root.node_id]
        is not previous_root_conditional
    )
    assert numpy.allclose(
        original_pass.conditional_by_node_id[left_internal.node_id],
        recomputed_pass.conditional_by_node_id[left_internal.node_id],
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.allclose(
        original_pass.conditional_by_node_id[tree.root.node_id],
        recomputed_pass.conditional_by_node_id[tree.root.node_id],
        rtol=0.0,
        atol=1e-12,
    )


def test_refresh_postorder_conditional_likelihoods_rejects_leaf_slice() -> None:
    tree = _balanced_four_tip_tree()
    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=2,
        leaf_likelihood=lambda _node: numpy.array([1.0, 0.0], dtype=float),
        transition_matrix_for_child=lambda _child: _transition_matrix(0.1),
    )

    with pytest.raises(ValueError, match="leaf node"):
        refresh_postorder_conditional_likelihoods(
            tree,
            pruning_pass,
            affected_node_ids=(tree.root.children[0].children[0].node_id,),
            transition_matrix_for_child=lambda _child: _transition_matrix(0.1),
        )


def _balanced_four_tip_tree() -> PhyloTree:
    return PhyloTree(
        TreeNode(
            children=[
                TreeNode(
                    branch_length=None,
                    children=[
                        TreeNode(name="A", branch_length=0.1),
                        TreeNode(name="B", branch_length=0.2),
                    ],
                ),
                TreeNode(
                    branch_length=None,
                    children=[
                        TreeNode(name="C", branch_length=0.3),
                        TreeNode(name="D", branch_length=0.4),
                    ],
                ),
            ]
        )
    )


def _leaf_vector(state: str) -> numpy.ndarray:
    if state == "0":
        return numpy.array([1.0, 0.0], dtype=float)
    return numpy.array([0.0, 1.0], dtype=float)


def _transition_matrix(branch_length: float) -> numpy.ndarray:
    stay = 1.0 - min(branch_length, 0.45)
    switch = 1.0 - stay
    return numpy.array(
        [[stay, switch], [switch / 2.0, 1.0 - (switch / 2.0)]],
        dtype=float,
    )
