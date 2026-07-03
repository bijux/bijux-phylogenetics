from __future__ import annotations

from dataclasses import dataclass
import math

import numpy

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


@dataclass(frozen=True, slots=True)
class ThreeTipPruningReference:
    """Closed-form three-tip pruning reference for one rooted A|(B,C) tree."""

    root_leaf_child: TreeNode
    internal_child: TreeNode
    internal_leaf_children: tuple[TreeNode, TreeNode]
    internal_raw: numpy.ndarray
    internal_normalized: numpy.ndarray
    internal_log_scale: float
    root_raw: numpy.ndarray
    root_normalized: numpy.ndarray
    root_log_scale: float
    log_likelihood: float


def build_three_tip_pruning_reference(
    *,
    tree: PhyloTree,
    root_prior: numpy.ndarray,
    leaf_vector_by_taxon: dict[str, numpy.ndarray],
    transition_matrix_by_node_id: dict[str, numpy.ndarray],
) -> ThreeTipPruningReference:
    """Build one exact three-tip pruning reference without using the runtime pass."""
    root_leaf_child, internal_child, internal_leaf_children = (
        _resolve_three_tip_fixture_shape(tree)
    )
    left_internal_leaf, right_internal_leaf = internal_leaf_children
    root_leaf_contribution = _transitioned_leaf_vector(
        root_leaf_child,
        leaf_vector_by_taxon=leaf_vector_by_taxon,
        transition_matrix_by_node_id=transition_matrix_by_node_id,
    )
    left_internal_contribution = _transitioned_leaf_vector(
        left_internal_leaf,
        leaf_vector_by_taxon=leaf_vector_by_taxon,
        transition_matrix_by_node_id=transition_matrix_by_node_id,
    )
    right_internal_contribution = _transitioned_leaf_vector(
        right_internal_leaf,
        leaf_vector_by_taxon=leaf_vector_by_taxon,
        transition_matrix_by_node_id=transition_matrix_by_node_id,
    )
    internal_raw = left_internal_contribution * right_internal_contribution
    internal_scale = float(internal_raw.sum())
    internal_normalized = internal_raw / internal_scale
    internal_log_scale = math.log(internal_scale)
    root_internal_contribution = transition_matrix_by_node_id[
        internal_child.node_id
    ] @ (internal_normalized)
    root_raw = root_leaf_contribution * root_internal_contribution
    root_scale = float(root_raw.sum())
    root_normalized = root_raw / root_scale
    root_log_scale = internal_log_scale + math.log(root_scale)
    log_likelihood = root_log_scale + math.log(
        float(numpy.dot(root_prior, root_normalized))
    )
    return ThreeTipPruningReference(
        root_leaf_child=root_leaf_child,
        internal_child=internal_child,
        internal_leaf_children=internal_leaf_children,
        internal_raw=internal_raw,
        internal_normalized=internal_normalized,
        internal_log_scale=internal_log_scale,
        root_raw=root_raw,
        root_normalized=root_normalized,
        root_log_scale=root_log_scale,
        log_likelihood=log_likelihood,
    )


def _resolve_three_tip_fixture_shape(
    tree: PhyloTree,
) -> tuple[TreeNode, TreeNode, tuple[TreeNode, TreeNode]]:
    root_children = tuple(tree.root.children)
    if len(root_children) != 2:
        raise ValueError("three-tip pruning reference requires a bifurcating root")
    leaf_children = tuple(child for child in root_children if child.is_leaf())
    internal_children = tuple(child for child in root_children if not child.is_leaf())
    if len(leaf_children) != 1 or len(internal_children) != 1:
        raise ValueError(
            "three-tip pruning reference requires a rooted A|(B,C) topology"
        )
    internal_child = internal_children[0]
    if len(internal_child.children) != 2 or any(
        not child.is_leaf() for child in internal_child.children
    ):
        raise ValueError(
            "three-tip pruning reference requires exactly two leaf descendants below the internal branch"
        )
    return leaf_children[0], internal_child, tuple(internal_child.children)


def _transitioned_leaf_vector(
    child: TreeNode,
    *,
    leaf_vector_by_taxon: dict[str, numpy.ndarray],
    transition_matrix_by_node_id: dict[str, numpy.ndarray],
) -> numpy.ndarray:
    if child.node_id is None:
        raise ValueError("tree node is missing a stable node_id")
    if child.name is None:
        raise ValueError("three-tip pruning reference requires named tips")
    return (
        transition_matrix_by_node_id[child.node_id] @ leaf_vector_by_taxon[child.name]
    )
