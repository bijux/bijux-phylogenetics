from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.metadata import load_taxon_table
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import MetadataJoinError
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class TreePruningReport:
    """Explicit record of tree pruning against a taxon table."""

    tree_path: Path
    keep_from_path: Path
    taxon_column: str
    original_tip_count: int
    kept_taxa: list[str]
    removed_taxa: list[str]


def _merge_branch_lengths(left: float | None, right: float | None) -> float | None:
    if left is None:
        return right
    if right is None:
        return left
    return left + right


def _collapse_unary(node: TreeNode) -> TreeNode:
    if len(node.children) != 1:
        return node
    only_child = _collapse_unary(node.children[0])
    return TreeNode(
        name=only_child.name,
        branch_length=_merge_branch_lengths(node.branch_length, only_child.branch_length),
        children=only_child.children,
    )


def _prune_node(node: TreeNode, keep_taxa: set[str]) -> TreeNode | None:
    if node.is_leaf():
        if node.name is None or node.name not in keep_taxa:
            return None
        return TreeNode(name=node.name, branch_length=node.branch_length, children=[])

    children = [child for child in (_prune_node(child, keep_taxa) for child in node.children) if child is not None]
    if not children:
        return None
    pruned_node = TreeNode(name=node.name, branch_length=node.branch_length, children=children)
    if len(children) == 1:
        return _collapse_unary(pruned_node)
    return pruned_node


def prune_tree_to_taxa(
    tree_path: Path,
    keep_from_path: Path,
    *,
    taxon_column: str | None = None,
) -> tuple[PhyloTree, TreePruningReport]:
    """Prune a tree to the taxa present in a metadata or traits table."""
    tree = load_tree(tree_path)
    keep_table = load_taxon_table(keep_from_path, taxon_column=taxon_column)
    keep_taxa = set(keep_table.taxa)
    retained_tips = sorted(name for name in tree.tip_names if name in keep_taxa)
    removed_tips = sorted(name for name in tree.tip_names if name not in keep_taxa)

    pruned_root = _prune_node(tree.root, keep_taxa)
    if pruned_root is None or not retained_tips:
        raise MetadataJoinError("no overlapping taxa remain after pruning request")

    pruned_root.branch_length = None
    pruned_tree = PhyloTree(root=pruned_root, source_format=tree.source_format)
    return pruned_tree, TreePruningReport(
        tree_path=tree_path,
        keep_from_path=keep_from_path,
        taxon_column=keep_table.taxon_column,
        original_tip_count=tree.tip_count,
        kept_taxa=retained_tips,
        removed_taxa=removed_tips,
    )


def write_pruned_taxa(path: Path, removed_taxa: list[str]) -> Path:
    """Write the removed taxa for an explicit pruning run."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["taxon"]
    lines.extend(removed_taxa)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
