from __future__ import annotations

from copy import deepcopy
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .models import TreeOrderingReport
from .transformation import (
    _compare_tree_topology,
    _copy_node_without_children,
    _descendant_taxa,
    _leaf_count,
    _summarize_transformation,
)


def _order_tree(node: TreeNode, *, strategy: str) -> TreeNode:
    if node.is_leaf():
        return _copy_node_without_children(node)

    ordered_children = [
        _order_tree(child, strategy=strategy) for child in node.children
    ]
    if strategy == "ladderize":
        ordered_children.sort(
            key=lambda child: (-_leaf_count(child), _descendant_taxa(child))
        )
    elif strategy == "alphabetical":
        ordered_children.sort(key=_descendant_taxa)
    else:
        raise ValueError(f"unsupported ordering strategy: {strategy}")

    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=ordered_children,
        metadata=deepcopy(node.metadata),
        edge_metadata=deepcopy(node.edge_metadata),
    )


def _rotate_named_node(
    node: TreeNode,
    *,
    clade_name: str,
) -> tuple[TreeNode, int]:
    rotated_children = []
    match_count = 1 if node.name == clade_name else 0
    for child in node.children:
        rotated_child, child_matches = _rotate_named_node(child, clade_name=clade_name)
        rotated_children.append(rotated_child)
        match_count += child_matches
    if node.name == clade_name:
        rotated_children = list(reversed(rotated_children))
    return (
        TreeNode(
            name=node.name,
            branch_length=node.branch_length,
            children=rotated_children,
            metadata=deepcopy(node.metadata),
            edge_metadata=deepcopy(node.edge_metadata),
        ),
        match_count,
    )


def _rotate_all_nodes(node: TreeNode) -> TreeNode:
    rotated_children = [_rotate_all_nodes(child) for child in node.children]
    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=list(reversed(rotated_children)),
        metadata=deepcopy(node.metadata),
        edge_metadata=deepcopy(node.edge_metadata),
    )


def ladderize_tree(tree_path: Path) -> tuple[PhyloTree, TreeOrderingReport]:
    """Ladderize a tree deterministically by descendant clade size."""
    tree = load_tree(tree_path)
    ladderized_tree = PhyloTree(
        root=_order_tree(tree.root, strategy="ladderize"),
        source_format=tree.source_format,
        rooted=tree.rooted,
    )
    comparison = _compare_tree_topology(tree, ladderized_tree)
    summary = _summarize_transformation(
        tree, ladderized_tree, transformation="ladderize-tree"
    )
    return ladderized_tree, TreeOrderingReport(
        tree_path=tree_path,
        strategy="ladderize",
        tip_order=ladderized_tree.tip_names,
        rooted_topology_preserved=comparison.topology_equal,
        unrooted_topology_preserved=comparison.same_unrooted_topology,
        summary=summary,
    )


def rotate_named_node(
    tree_path: Path, *, clade_name: str
) -> tuple[PhyloTree, TreeOrderingReport]:
    """Reverse the child order at one named internal node."""
    tree = load_tree(tree_path)
    rotated_root, match_count = _rotate_named_node(tree.root, clade_name=clade_name)
    if match_count == 0:
        raise ValueError(f"clade '{clade_name}' was not found in {tree_path}")
    if match_count > 1:
        raise ValueError(f"clade '{clade_name}' is ambiguous in {tree_path}")
    rotated_tree = PhyloTree(
        root=rotated_root,
        source_format=tree.source_format,
        rooted=tree.rooted,
    )
    comparison = _compare_tree_topology(tree, rotated_tree)
    summary = _summarize_transformation(
        tree,
        rotated_tree,
        transformation="rotate-named-node",
        extra_changed_nodes=[clade_name],
    )
    return rotated_tree, TreeOrderingReport(
        tree_path=tree_path,
        strategy=f"rotate:{clade_name}",
        tip_order=rotated_tree.tip_names,
        rooted_topology_preserved=comparison.topology_equal,
        unrooted_topology_preserved=comparison.same_unrooted_topology,
        summary=summary,
    )


def rotate_all_internal_nodes(tree_path: Path) -> tuple[PhyloTree, TreeOrderingReport]:
    """Reverse the child order at every internal node."""
    tree = load_tree(tree_path)
    rotated_tree = PhyloTree(
        root=_rotate_all_nodes(tree.root),
        source_format=tree.source_format,
        rooted=tree.rooted,
    )
    comparison = _compare_tree_topology(tree, rotated_tree)
    summary = _summarize_transformation(
        tree,
        rotated_tree,
        transformation="rotate-all-internal-nodes",
    )
    return rotated_tree, TreeOrderingReport(
        tree_path=tree_path,
        strategy="rotate-all",
        tip_order=rotated_tree.tip_names,
        rooted_topology_preserved=comparison.topology_equal,
        unrooted_topology_preserved=comparison.same_unrooted_topology,
        summary=summary,
    )


def sort_tree_tips_alphabetically(
    tree_path: Path,
) -> tuple[PhyloTree, TreeOrderingReport]:
    """Sort tree children recursively by alphabetical descendant tip order."""
    tree = load_tree(tree_path)
    sorted_tree = PhyloTree(
        root=_order_tree(tree.root, strategy="alphabetical"),
        source_format=tree.source_format,
        rooted=tree.rooted,
    )
    comparison = _compare_tree_topology(tree, sorted_tree)
    summary = _summarize_transformation(
        tree, sorted_tree, transformation="sort-tree-tips"
    )
    return sorted_tree, TreeOrderingReport(
        tree_path=tree_path,
        strategy="alphabetical",
        tip_order=sorted_tree.tip_names,
        rooted_topology_preserved=comparison.topology_equal,
        unrooted_topology_preserved=comparison.same_unrooted_topology,
        summary=summary,
    )
