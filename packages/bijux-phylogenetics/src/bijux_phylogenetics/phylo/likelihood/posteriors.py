from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .pruning import postorder_conditional_likelihoods


@dataclass(slots=True)
class FiniteStateMarginalPosteriorPass:
    """Marginal posterior vectors for one finite-state rooted tree likelihood."""

    posterior_by_node_id: dict[str, numpy.ndarray]
    outside_by_node_id: dict[str, numpy.ndarray]
    state_count: int

    def posterior_for_node(self, node: TreeNode) -> numpy.ndarray:
        if node.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        return self.posterior_by_node_id[node.node_id]

    def outside_for_node(self, node: TreeNode) -> numpy.ndarray:
        if node.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        return self.outside_by_node_id[node.node_id]


def compute_marginal_state_posteriors(
    tree: PhyloTree,
    *,
    state_count: int,
    leaf_likelihood: Callable[[TreeNode], numpy.ndarray],
    transition_matrix_for_child: Callable[[TreeNode], numpy.ndarray],
    root_prior: numpy.ndarray,
) -> FiniteStateMarginalPosteriorPass:
    """Compute node-wise marginal state posteriors for one finite-state site."""
    if state_count < 1:
        raise ValueError("state_count must be positive")

    normalized_root_prior = _normalize_probability_vector(
        root_prior,
        state_count=state_count,
        context="root prior",
    )
    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=state_count,
        leaf_likelihood=leaf_likelihood,
        transition_matrix_for_child=transition_matrix_for_child,
        normalize_internal_conditionals=True,
    )
    posterior_by_node_id: dict[str, numpy.ndarray] = {}
    outside_by_node_id: dict[str, numpy.ndarray] = {}
    root_id = _node_id(tree.root)
    outside_by_node_id[root_id] = normalized_root_prior
    posterior_by_node_id[root_id] = _normalize_probability_vector(
        normalized_root_prior * pruning_pass.conditional_for_node(tree.root),
        state_count=state_count,
        context=f"node '{root_id}' posterior",
    )

    for node in tree.iter_nodes(order="preorder"):
        node_id = _node_id(node)
        parent_outside = outside_by_node_id[node_id]
        if node.is_leaf():
            continue
        child_messages = {
            _node_id(child): _transition_matrix(
                transition_matrix_for_child(child),
                state_count=state_count,
                context=f"branch '{_node_id(child)}'",
            )
            @ pruning_pass.conditional_for_node(child)
            for child in node.children
        }
        for child in node.children:
            child_id = _node_id(child)
            parent_context = parent_outside.copy()
            for sibling in node.children:
                sibling_id = _node_id(sibling)
                if sibling_id == child_id:
                    continue
                parent_context *= child_messages[sibling_id]
            branch_transition = _transition_matrix(
                transition_matrix_for_child(child),
                state_count=state_count,
                context=f"branch '{child_id}'",
            )
            child_outside = _normalize_probability_vector(
                branch_transition.T @ parent_context,
                state_count=state_count,
                context=f"node '{child_id}' outside probabilities",
            )
            outside_by_node_id[child_id] = child_outside
            posterior_by_node_id[child_id] = _normalize_probability_vector(
                child_outside * pruning_pass.conditional_for_node(child),
                state_count=state_count,
                context=f"node '{child_id}' posterior",
            )

    return FiniteStateMarginalPosteriorPass(
        posterior_by_node_id=posterior_by_node_id,
        outside_by_node_id=outside_by_node_id,
        state_count=state_count,
    )


def _node_id(node: TreeNode) -> str:
    if node.node_id is None:
        raise ValueError("tree node is missing a stable node_id")
    return node.node_id


def _normalize_probability_vector(
    vector: numpy.ndarray,
    *,
    state_count: int,
    context: str,
) -> numpy.ndarray:
    candidate = numpy.asarray(vector, dtype=float)
    if candidate.shape != (state_count,):
        raise ValueError(
            f"{context} must be a one-dimensional vector with {state_count} states"
        )
    if not numpy.all(numpy.isfinite(candidate)):
        raise ValueError(f"{context} must contain only finite values")
    if numpy.any(candidate < 0.0):
        raise ValueError(f"{context} must not contain negative probabilities")
    total = float(candidate.sum())
    if total <= 0.0:
        raise ValueError(f"{context} must sum to a positive value")
    return candidate / total


def _transition_matrix(
    matrix: numpy.ndarray,
    *,
    state_count: int,
    context: str,
) -> numpy.ndarray:
    candidate = numpy.asarray(matrix, dtype=float)
    if candidate.shape != (state_count, state_count):
        raise ValueError(
            f"{context} must be a {state_count}x{state_count} transition matrix"
        )
    if not numpy.all(numpy.isfinite(candidate)):
        raise ValueError(f"{context} must contain only finite values")
    if numpy.any(candidate < 0.0):
        raise ValueError(f"{context} must not contain negative probabilities")
    return candidate.copy()
