from __future__ import annotations

from pathlib import Path
import re

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import InvalidBranchLengthError
from bijux_phylogenetics.io.biopython import load_biophylo, loads_biophylo

_BRANCH_LENGTH_PATTERN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$")
_UNQUOTED_LABEL_PATTERN = re.compile(r"^[0-9A-Za-z._-]+$")


def _validate_branch_lengths(text: str) -> None:
    index = 0
    while True:
        index = text.find(":", index)
        if index == -1:
            return
        cursor = index + 1
        while cursor < len(text) and text[cursor].isspace():
            cursor += 1
        start = cursor
        while cursor < len(text) and text[cursor] not in ",();":
            cursor += 1
        raw_value = text[start:cursor].strip()
        if not raw_value:
            raise InvalidBranchLengthError(
                f"missing branch length in Newick string at position {start}"
            )
        if not _BRANCH_LENGTH_PATTERN.fullmatch(raw_value):
            raise InvalidBranchLengthError(
                f"invalid branch length '{raw_value}' in Newick string"
            )
        index = cursor


def loads_newick(text: str):
    """Parse a Newick string into a minimal tree object."""
    _validate_branch_lengths(text)
    return loads_biophylo(text, source_format="newick")


def load_newick(path: Path):
    """Load a Newick tree from disk."""
    return load_biophylo(path, source_format="newick")


def _sort_key(node: TreeNode) -> tuple[str, int]:
    tip_names = sorted(name for name in PhyloTree(root=node).tip_names if name)
    return (tip_names[0] if tip_names else "", len(tip_names))


def _format_branch_length(branch_length: float | None) -> str:
    if branch_length is None:
        return ""
    return f":{format(branch_length, '.15g')}"


def _format_label(label: str | None) -> str:
    if label is None:
        return ""
    if _UNQUOTED_LABEL_PATTERN.fullmatch(label):
        return label
    escaped = label.replace("'", "''")
    return f"'{escaped}'"


def _serialize_node(node: TreeNode) -> str:
    if node.children:
        ordered_children = ",".join(
            _serialize_node(child) for child in sorted(node.children, key=_sort_key)
        )
        label = _format_label(node.name)
        return f"({ordered_children}){label}{_format_branch_length(node.branch_length)}"
    return f"{_format_label(node.name)}{_format_branch_length(node.branch_length)}"


def dumps_newick(tree: PhyloTree) -> str:
    """Serialize a local tree model into canonical Newick."""
    if tree.root.children:
        ordered_children = ",".join(
            _serialize_node(child)
            for child in sorted(tree.root.children, key=_sort_key)
        )
        root_label = _format_label(tree.root.name)
        return f"({ordered_children}){root_label};"
    return f"{_serialize_node(tree.root)};"


def write_newick(path: Path, tree: PhyloTree) -> Path:
    """Write a canonical Newick tree to disk."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"{dumps_newick(tree)}\n", encoding="utf-8")
    return path
