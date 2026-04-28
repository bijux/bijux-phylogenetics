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
