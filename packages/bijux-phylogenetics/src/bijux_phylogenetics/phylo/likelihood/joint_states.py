from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


@dataclass(slots=True)
class FiniteStateJointAssignmentPass:
    """One max-product joint state assignment across one rooted finite-state tree."""

    assigned_state_index_by_node_id: dict[str, int]
    conditioned_subtree_score_by_node_id: dict[str, numpy.ndarray]
    max_joint_score: float
    root_state_index: int
    state_count: int

    def assigned_state_index_for_node(self, node: TreeNode) -> int:
        if node.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        return self.assigned_state_index_by_node_id[node.node_id]

    def conditioned_subtree_score_for_node(self, node: TreeNode) -> numpy.ndarray:
        if node.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        return self.conditioned_subtree_score_by_node_id[node.node_id]


def compute_joint_state_assignment(
    tree: PhyloTree,
    *,
    state_count: int,
    leaf_likelihood: Callable[[TreeNode], numpy.ndarray],
    transition_matrix_for_child: Callable[[TreeNode], numpy.ndarray],
    root_prior: numpy.ndarray,
) -> FiniteStateJointAssignmentPass:
    """Compute one globally optimal rooted state assignment by max-product recursion."""
    if state_count < 1:
        raise ValueError("state_count must be positive")

    validated_root_prior = _validated_nonnegative_vector(
        root_prior,
        state_count=state_count,
        context="root prior",
    )
    score_by_node_id: dict[str, numpy.ndarray] = {}
    argmax_state_by_child_id: dict[str, numpy.ndarray] = {}

    for node in tree.iter_nodes(order="postorder"):
        node_id = _node_id(node)
        if node.is_leaf():
            score_by_node_id[node_id] = _validated_nonnegative_vector(
                leaf_likelihood(node),
                state_count=state_count,
                context=f"leaf '{node.name or node_id}'",
            )
            continue

        conditioned_scores = numpy.ones(state_count, dtype=float)
        for child in node.children:
            child_id = _node_id(child)
            child_scores = score_by_node_id[child_id]
            transition = _validated_transition_matrix(
                transition_matrix_for_child(child),
                state_count=state_count,
                context=f"branch '{child_id}'",
            )
            best_child_scores = numpy.zeros(state_count, dtype=float)
            best_child_states = numpy.zeros(state_count, dtype=int)
            for parent_state_index in range(state_count):
                weighted_scores = transition[parent_state_index, :] * child_scores
                best_child_state_index = _stable_argmax(weighted_scores)
                best_child_scores[parent_state_index] = float(
                    weighted_scores[best_child_state_index]
                )
                best_child_states[parent_state_index] = best_child_state_index
            conditioned_scores *= best_child_scores
            argmax_state_by_child_id[child_id] = best_child_states
        score_by_node_id[node_id] = conditioned_scores

    root_scores = validated_root_prior * score_by_node_id[_node_id(tree.root)]
    root_state_index = _stable_argmax(root_scores)
    assigned_state_index_by_node_id = {_node_id(tree.root): root_state_index}
    _traceback_joint_state_assignment(
        tree.root,
        root_state_index=root_state_index,
        assigned_state_index_by_node_id=assigned_state_index_by_node_id,
        argmax_state_by_child_id=argmax_state_by_child_id,
    )

    return FiniteStateJointAssignmentPass(
        assigned_state_index_by_node_id=assigned_state_index_by_node_id,
        conditioned_subtree_score_by_node_id=score_by_node_id,
        max_joint_score=float(root_scores[root_state_index]),
        root_state_index=root_state_index,
        state_count=state_count,
    )


def _traceback_joint_state_assignment(
    node: TreeNode,
    *,
    root_state_index: int,
    assigned_state_index_by_node_id: dict[str, int],
    argmax_state_by_child_id: dict[str, numpy.ndarray],
) -> None:
    if node.is_leaf():
        return
    for child in node.children:
        child_id = _node_id(child)
        child_state_index = int(argmax_state_by_child_id[child_id][root_state_index])
        assigned_state_index_by_node_id[child_id] = child_state_index
        if not child.is_leaf():
            _traceback_joint_state_assignment(
                child,
                root_state_index=child_state_index,
                assigned_state_index_by_node_id=assigned_state_index_by_node_id,
                argmax_state_by_child_id=argmax_state_by_child_id,
            )


def _stable_argmax(values: numpy.ndarray) -> int:
    best_index = 0
    best_value = float(values[0])
    for index in range(1, len(values)):
        candidate_value = float(values[index])
        if candidate_value > best_value:
            best_value = candidate_value
            best_index = index
    return best_index


def _node_id(node: TreeNode) -> str:
    if node.node_id is None:
        raise ValueError("tree node is missing a stable node_id")
    return node.node_id


def _validated_nonnegative_vector(
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
        raise ValueError(f"{context} must not contain negative values")
    if float(candidate.sum()) <= 0.0:
        raise ValueError(f"{context} must contain at least one positive value")
    return candidate.copy()


def _validated_transition_matrix(
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
