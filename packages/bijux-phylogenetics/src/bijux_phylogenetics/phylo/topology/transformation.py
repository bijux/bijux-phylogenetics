from __future__ import annotations

from copy import deepcopy
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.clades import (
    informative_rooted_clades,
    informative_unrooted_splits,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .models import BranchCollapseReport, TreeTransformationSummary, _TopologyComparison


def _clone_node(node: TreeNode) -> TreeNode:
    return node.copy()


def _copy_node_without_children(node: TreeNode) -> TreeNode:
    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=[],
        metadata=deepcopy(node.metadata),
        edge_metadata=deepcopy(node.edge_metadata),
    )


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return sorted(taxa)


def _join_taxa(taxa: list[str]) -> str:
    return ",".join(taxa)


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(round(value, 15))


def _leaf_count(node: TreeNode) -> int:
    if node.is_leaf():
        return 1
    return sum(_leaf_count(child) for child in node.children)


def _combine_branch_lengths(base: float | None, extra: float | None) -> float | None:
    if base is None and extra is None:
        return None
    if base is None:
        return extra
    if extra is None:
        return base
    return base + extra


def _node_signature(node: TreeNode) -> str:
    taxa = _descendant_taxa(node)
    if taxa:
        return "|".join(taxa)
    return node.name or "<unnamed>"


def _branch_length_affecting_nodes(
    original: TreeNode, transformed: TreeNode
) -> list[str]:
    original_nodes = {
        _node_signature(node): node.branch_length
        for node in original.iter_nodes()
        if node is not original
    }
    transformed_nodes = {
        _node_signature(node): node.branch_length
        for node in transformed.iter_nodes()
        if node is not transformed
    }
    affected = {
        signature
        for signature in set(original_nodes) | set(transformed_nodes)
        if original_nodes.get(signature) != transformed_nodes.get(signature)
    }
    return sorted(affected)


def _changed_node_signatures(original: TreeNode, transformed: TreeNode) -> list[str]:
    original_nodes = {_node_signature(node) for node in original.iter_nodes()}
    transformed_nodes = {_node_signature(node) for node in transformed.iter_nodes()}
    return sorted(original_nodes.symmetric_difference(transformed_nodes))


def _summarize_transformation(
    original: PhyloTree,
    transformed: PhyloTree,
    *,
    transformation: str,
    extra_changed_nodes: list[str] | None = None,
) -> TreeTransformationSummary:
    original_taxa = sorted(original.tip_names)
    transformed_taxa = sorted(transformed.tip_names)
    original_taxa_set = set(original_taxa)
    transformed_taxa_set = set(transformed_taxa)
    nodes_changed = _changed_node_signatures(original.root, transformed.root)
    if extra_changed_nodes:
        nodes_changed = sorted(set(nodes_changed) | set(extra_changed_nodes))
    branch_lengths_affected = _branch_length_affecting_nodes(
        original.root, transformed.root
    )
    original_total = round(original.total_branch_length(), 15)
    transformed_total = round(transformed.total_branch_length(), 15)
    return TreeTransformationSummary(
        transformation=transformation,
        original_tip_count=original.tip_count,
        transformed_tip_count=transformed.tip_count,
        retained_taxa=sorted(original_taxa_set & transformed_taxa_set),
        removed_taxa=sorted(original_taxa_set - transformed_taxa_set),
        added_taxa=sorted(transformed_taxa_set - original_taxa_set),
        original_internal_node_count=original.internal_node_count,
        transformed_internal_node_count=transformed.internal_node_count,
        nodes_changed=nodes_changed,
        original_total_branch_length=original_total,
        transformed_total_branch_length=transformed_total,
        branch_length_delta=round(transformed_total - original_total, 15),
        branch_lengths_affected=branch_lengths_affected,
    )


def _compare_tree_topology(
    original: PhyloTree, transformed: PhyloTree
) -> _TopologyComparison:
    shared_taxa = set(original.tip_names) & set(transformed.tip_names)
    if len(shared_taxa) < 2:
        return _TopologyComparison(
            topology_equal=original.tip_names == transformed.tip_names,
            same_unrooted_topology=True,
        )
    left_clades = informative_rooted_clades(original, shared_taxa)
    right_clades = informative_rooted_clades(transformed, shared_taxa)
    return _TopologyComparison(
        topology_equal=left_clades == right_clades,
        same_unrooted_topology=informative_unrooted_splits(original, shared_taxa)
        == informative_unrooted_splits(transformed, shared_taxa),
    )


def _collapse_short_internal_branches(
    node: TreeNode,
    *,
    threshold: float,
    collapsed_clades: list[str],
) -> TreeNode:
    if node.is_leaf():
        return _copy_node_without_children(node)

    rewritten_children: list[TreeNode] = []
    for child in node.children:
        rewritten_child = _collapse_short_internal_branches(
            child,
            threshold=threshold,
            collapsed_clades=collapsed_clades,
        )
        should_collapse = (
            not rewritten_child.is_leaf()
            and rewritten_child.branch_length is not None
            and rewritten_child.branch_length < threshold
        )
        if should_collapse:
            collapsed_clades.append(_node_signature(rewritten_child))
            rewritten_children.extend(
                TreeNode(
                    name=grandchild.name,
                    branch_length=grandchild.branch_length,
                    children=[
                        _clone_node(child_node) for child_node in grandchild.children
                    ],
                    metadata=deepcopy(grandchild.metadata),
                    edge_metadata=deepcopy(grandchild.edge_metadata),
                )
                for grandchild in rewritten_child.children
            )
            continue
        rewritten_children.append(rewritten_child)

    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=rewritten_children,
        metadata=deepcopy(node.metadata),
        edge_metadata=deepcopy(node.edge_metadata),
    )


def collapse_branches_below_length(
    tree_path: Path,
    *,
    threshold: float,
) -> tuple[PhyloTree, BranchCollapseReport]:
    """Collapse short internal branches into parent-level polytomies."""
    if threshold < 0:
        raise ValueError(f"threshold must be non-negative, got {threshold}")

    tree = load_tree(tree_path)
    collapsed_clades: list[str] = []
    collapsed_root = _collapse_short_internal_branches(
        tree.root,
        threshold=threshold,
        collapsed_clades=collapsed_clades,
    )
    collapsed_root.branch_length = None
    collapsed_tree = PhyloTree(
        root=collapsed_root, source_format=tree.source_format, rooted=tree.rooted
    )
    comparison = _compare_tree_topology(tree, collapsed_tree)
    summary = _summarize_transformation(
        tree,
        collapsed_tree,
        transformation="collapse-short-branches",
        extra_changed_nodes=collapsed_clades,
    )
    return collapsed_tree, BranchCollapseReport(
        tree_path=tree_path,
        threshold=threshold,
        collapsed_clades=sorted(collapsed_clades),
        topology_preserved=comparison.same_unrooted_topology,
        summary=summary,
    )
