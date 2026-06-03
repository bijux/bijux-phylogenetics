from __future__ import annotations

from io import StringIO
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import TreeParseError


def _convert_clade(clade: Clade) -> TreeNode:
    label = clade.name
    if label is None and clade.confidence is not None:
        label = format(float(clade.confidence), ".15g")
    metadata: dict[str, object] = {}
    if clade.confidence is not None:
        metadata["confidence"] = float(clade.confidence)
    comment = getattr(clade, "comment", None)
    if comment is not None:
        metadata["comment"] = comment
    return TreeNode(
        name=label,
        branch_length=clade.branch_length,
        children=[_convert_clade(child) for child in clade.clades],
        metadata=metadata,
    )


def _convert_tree_node(node: TreeNode) -> Clade:
    return Clade(
        branch_length=node.branch_length,
        name=node.name,
        confidence=(
            float(node.metadata["confidence"])
            if "confidence" in node.metadata
            else None
        ),
        clades=[_convert_tree_node(child) for child in node.children],
    )


def tree_from_biophylo(tree: Tree, *, source_format: str) -> PhyloTree:
    """Convert a Bio.Phylo tree into the local tree model."""
    return PhyloTree(
        root=_convert_clade(tree.root), source_format=source_format, rooted=tree.rooted
    )


def tree_to_biophylo(tree: PhyloTree) -> Tree:
    """Convert the local tree model into a Bio.Phylo tree."""
    return Tree(root=_convert_tree_node(tree.root), rooted=bool(tree.rooted))


def loads_biophylo(text: str, *, source_format: str) -> PhyloTree:
    """Parse tree text with Bio.Phylo and convert it to the local tree model."""
    try:
        tree = Phylo.read(StringIO(text), source_format)
    except Exception as error:  # pragma: no cover - biopython exception surface varies
        raise TreeParseError(str(error)) from error
    return tree_from_biophylo(tree, source_format=source_format)


def load_biophylo(path: Path, *, source_format: str) -> PhyloTree:
    """Read a tree file with Bio.Phylo and convert it to the local tree model."""
    if not path.exists():
        raise FileNotFoundError(f"tree file not found: {path}")
    try:
        tree = Phylo.read(path, source_format)
    except Exception as error:  # pragma: no cover - biopython exception surface varies
        raise TreeParseError(str(error)) from error
    return tree_from_biophylo(tree, source_format=source_format)
