from __future__ import annotations

import math
from pathlib import Path
import re

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    TreeParseError,
    UnnamedTipError,
)

_BRANCH_LENGTH_PATTERN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$")
_UNQUOTED_LABEL_PATTERN = re.compile(r"^[0-9A-Za-z._-]+$")
_NUMERIC_LABEL_PATTERN = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$")


class _NewickParser:
    """Recursive-descent Newick parser for the native `PhyloTree` model."""

    def __init__(self, text: str) -> None:
        self.text = text
        self.index = 0

    def parse_tree(self) -> PhyloTree:
        self._skip_whitespace()
        root = self._parse_subtree()
        self._skip_whitespace()
        self._expect(";")
        return PhyloTree(root=root, source_format="newick", rooted=False)

    def parse_tree_set(self) -> list[PhyloTree]:
        trees: list[PhyloTree] = []
        while True:
            self._skip_whitespace()
            if self._at_end():
                break
            trees.append(self.parse_tree())
        if not trees:
            raise TreeParseError("tree set contains no trees")
        return trees

    def _parse_subtree(self) -> TreeNode:
        self._skip_whitespace()
        if self._peek() == "(":
            return self._parse_internal_node()
        return self._parse_leaf()

    def _parse_internal_node(self) -> TreeNode:
        self._expect("(")
        children: list[TreeNode] = []
        while True:
            self._skip_whitespace()
            if self._peek() == ")":
                if not children:
                    raise self._parse_error("empty internal node is not valid Newick")
                self._advance()
                break
            children.append(self._parse_subtree())
            self._skip_whitespace()
            delimiter = self._peek()
            if delimiter == ",":
                self._advance()
                continue
            if delimiter == ")":
                self._advance()
                break
            if delimiter is None:
                raise self._parse_error(
                    "unexpected end of Newick string while parsing child list"
                )
            raise self._parse_error(
                f"unexpected token '{delimiter}' while parsing child list"
            )
        label = self._parse_optional_label()
        branch_length = self._parse_optional_branch_length()
        return self._build_node(
            children=children,
            label=label,
            branch_length=branch_length,
        )

    def _parse_leaf(self) -> TreeNode:
        label = self._parse_optional_label()
        branch_length = self._parse_optional_branch_length()
        if label is None and branch_length is None:
            raise self._parse_error("terminal node is missing a label")
        return TreeNode(name=label, branch_length=branch_length)

    def _parse_optional_label(self) -> str | None:
        self._skip_whitespace()
        token = self._peek()
        if token is None or token in ",():;":  # nosec B105
            return None
        if token == "'":  # nosec B105
            return self._parse_quoted_label()
        return self._parse_unquoted_label()

    def _parse_quoted_label(self) -> str:
        self._expect("'")
        chunks: list[str] = []
        while True:
            token = self._peek()
            if token is None:
                raise self._parse_error("unterminated quoted label")
            if token == "'":  # nosec B105
                self._advance()
                if self._peek() == "'":
                    self._advance()
                    chunks.append("'")
                    continue
                break
            chunks.append(self._advance())
        return "".join(chunks)

    def _parse_unquoted_label(self) -> str:
        start = self.index
        while True:
            token = self._peek()
            if token is None or token in ",():;":  # nosec B105
                break
            self._advance()
        label = self.text[start : self.index].strip()
        if not label:
            raise self._parse_error("node label is empty")
        return label

    def _parse_optional_branch_length(self) -> float | None:
        self._skip_whitespace()
        if self._peek() != ":":
            return None
        self._advance()
        self._skip_whitespace()
        start = self.index
        while True:
            token = self._peek()
            if token is None or token in ",();":  # nosec B105
                break
            self._advance()
        raw_value = self.text[start : self.index].strip()
        if not raw_value:
            raise self._branch_length_error("missing branch length", start)
        if not _BRANCH_LENGTH_PATTERN.fullmatch(raw_value):
            raise self._branch_length_error(
                f"invalid branch length '{raw_value}'",
                start,
            )
        return float(raw_value)

    def _build_node(
        self,
        *,
        children: list[TreeNode],
        label: str | None,
        branch_length: float | None,
    ) -> TreeNode:
        metadata: dict[str, object] = {}
        if label is not None and _NUMERIC_LABEL_PATTERN.fullmatch(label):
            metadata["confidence"] = float(label)
        return TreeNode(
            name=label,
            branch_length=branch_length,
            children=children,
            metadata=metadata,
        )

    def _expect(self, token: str) -> None:
        self._skip_whitespace()
        observed = self._peek()
        if observed != token:
            if observed is None:
                raise self._parse_error(
                    f"expected '{token}' before end of Newick string"
                )
            raise self._parse_error(f"expected '{token}' but found '{observed}'")
        self._advance()

    def _skip_whitespace(self) -> None:
        while not self._at_end() and self.text[self.index].isspace():
            self.index += 1

    def _peek(self) -> str | None:
        if self._at_end():
            return None
        return self.text[self.index]

    def _advance(self) -> str:
        token = self.text[self.index]
        self.index += 1
        return token

    def _at_end(self) -> bool:
        return self.index >= len(self.text)

    def _line_column(self, position: int) -> tuple[int, int]:
        line = self.text.count("\n", 0, position) + 1
        line_start = self.text.rfind("\n", 0, position) + 1
        column = (position - line_start) + 1
        return line, column

    def _parse_error(self, message: str, position: int | None = None) -> TreeParseError:
        error_position = self.index if position is None else position
        line, column = self._line_column(error_position)
        return TreeParseError(
            f"{message} at line {line}, column {column}",
            details={
                "position": error_position,
                "line": line,
                "column": column,
            },
        )

    def _branch_length_error(
        self,
        message: str,
        position: int,
    ) -> InvalidBranchLengthError:
        line, column = self._line_column(position)
        return InvalidBranchLengthError(
            f"{message} at line {line}, column {column}",
            details={
                "position": position,
                "line": line,
                "column": column,
            },
        )


def loads_newick(text: str):
    """Parse a Newick string into a minimal tree object."""
    parser = _NewickParser(text)
    tree = parser.parse_tree()
    parser._skip_whitespace()
    if not parser._at_end():
        raise parser._parse_error("unexpected trailing content after Newick record")
    return tree


def load_newick(path: Path):
    """Load a Newick tree from disk."""
    if not path.exists():
        raise FileNotFoundError(f"tree file not found: {path}")
    return loads_newick(path.read_text(encoding="utf-8"))


def loads_newick_tree_set(text: str) -> list[PhyloTree]:
    """Load one Newick tree per record from a text buffer."""
    return _NewickParser(text).parse_tree_set()


def load_newick_tree_set(path: Path) -> list[PhyloTree]:
    """Load one Newick tree per record from disk."""
    if not path.exists():
        raise FileNotFoundError(f"tree-set file not found: {path}")
    return loads_newick_tree_set(path.read_text(encoding="utf-8"))


def iter_newick_tree_records(text: str):
    """Yield one normalized Newick record per parsed statement."""
    record_index = 0
    buffer = ""
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        buffer = f"{buffer} {line}".strip() if buffer else line
        while ";" in buffer:
            statement, buffer = buffer.split(";", 1)
            normalized = statement.strip()
            if not normalized:
                continue
            record_index += 1
            yield record_index, f"{normalized};"
    remainder = buffer.strip()
    if remainder:
        record_index += 1
        yield record_index, remainder


def iter_newick_tree_records_from_path(path: Path):
    """Yield one normalized Newick record per parsed statement from disk."""
    if not path.exists():
        raise FileNotFoundError(f"tree-set file not found: {path}")
    yield from iter_newick_tree_records(path.read_text(encoding="utf-8"))


def quote_newick_label(label: str) -> str:
    """Quote one tip or node label for safe Newick reuse."""
    return _format_label(label)


def _sort_key(node: TreeNode) -> tuple[str, int]:
    tip_names = sorted(name for name in node.descendant_taxa if name)
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


def _validate_tree_for_newick(node: TreeNode) -> None:
    if node.is_leaf() and node.name is None:
        raise UnnamedTipError(
            "tree contains unnamed tips and cannot be written as Newick"
        )
    if node.branch_length is not None and not math.isfinite(node.branch_length):
        raise InvalidBranchLengthError(
            f"invalid branch length {node.branch_length!r} in tree node"
        )
    for child in node.children:
        _validate_tree_for_newick(child)


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
    _validate_tree_for_newick(tree.root)
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


def dumps_newick_tree_set(trees: list[PhyloTree]) -> str:
    """Serialize a list of trees as one canonical Newick record per line."""
    if not trees:
        raise TreeParseError(
            "tree set contains no trees and cannot be written as Newick"
        )
    return "".join(f"{dumps_newick(tree)}\n" for tree in trees)


def write_newick_tree_set(path: Path, trees: list[PhyloTree]) -> Path:
    """Write one canonical Newick tree per line."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dumps_newick_tree_set(trees), encoding="utf-8")
    return path
