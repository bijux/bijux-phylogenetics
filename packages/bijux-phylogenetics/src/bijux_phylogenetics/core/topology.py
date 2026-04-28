from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class CladeExtractionReport:
    """Explicit record of subtree extraction for a named clade."""

    tree_path: Path
    clade_name: str
    tip_count: int
    taxa: list[str]


@dataclass(slots=True)
class BranchCollapseReport:
    """Explicit record of internal branches collapsed by a length threshold."""

    tree_path: Path
    threshold: float
    collapsed_clades: list[str]


def _clone_node(node: TreeNode) -> TreeNode:
    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=[_clone_node(child) for child in node.children],
    )


def _find_named_nodes(node: TreeNode, *, clade_name: str) -> list[TreeNode]:
    matches: list[TreeNode] = []
    if node.name == clade_name:
        matches.append(node)
    for child in node.children:
        matches.extend(_find_named_nodes(child, clade_name=clade_name))
    return matches


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return sorted(taxa)


def _node_signature(node: TreeNode) -> str:
    taxa = _descendant_taxa(node)
    if taxa:
        return "|".join(taxa)
    return node.name or "<unnamed>"


def _collapse_short_internal_branches(
    node: TreeNode,
    *,
    threshold: float,
    collapsed_clades: list[str],
) -> TreeNode:
    if node.is_leaf():
        return TreeNode(name=node.name, branch_length=node.branch_length, children=[])

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
                    children=[_clone_node(child_node) for child_node in grandchild.children],
                )
                for grandchild in rewritten_child.children
            )
            continue
        rewritten_children.append(rewritten_child)

    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=rewritten_children,
    )


def extract_named_clade(
    tree_path: Path,
    *,
    clade_name: str,
) -> tuple[PhyloTree, CladeExtractionReport]:
    """Extract a named clade as a standalone subtree."""
    tree = load_tree(tree_path)
    matches = _find_named_nodes(tree.root, clade_name=clade_name)
    if not matches:
        raise ValueError(f"clade '{clade_name}' was not found in {tree_path}")
    if len(matches) > 1:
        raise ValueError(f"clade '{clade_name}' is ambiguous in {tree_path}")

    subtree_root = _clone_node(matches[0])
    subtree_root.branch_length = None
    subtree = PhyloTree(root=subtree_root, source_format=tree.source_format)
    return subtree, CladeExtractionReport(
        tree_path=tree_path,
        clade_name=clade_name,
        tip_count=subtree.tip_count,
        taxa=sorted(subtree.tip_names),
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
    collapsed_tree = PhyloTree(root=collapsed_root, source_format=tree.source_format)
    return collapsed_tree, BranchCollapseReport(
        tree_path=tree_path,
        threshold=threshold,
        collapsed_clades=sorted(collapsed_clades),
    )
