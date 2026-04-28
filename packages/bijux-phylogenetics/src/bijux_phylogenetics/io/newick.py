from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode


class _NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.length = len(self.text)
        self.index = 0

    def parse(self) -> PhyloTree:
        node = self._parse_subtree()
        self._skip_whitespace()
        if self._peek() == ";":
            self.index += 1
        self._skip_whitespace()
        if self.index != self.length:
            raise ValueError(f"unexpected trailing content in Newick string at position {self.index}")
        return PhyloTree(root=node, source_format="newick")

    def _parse_subtree(self) -> TreeNode:
        self._skip_whitespace()
        if self._peek() == "(":
            self.index += 1
            children: list[TreeNode] = []
            while True:
                children.append(self._parse_subtree())
                self._skip_whitespace()
                token = self._peek()
                if token == ",":
                    self.index += 1
                    continue
                if token == ")":
                    self.index += 1
                    break
                raise ValueError(f"expected ',' or ')' in Newick string at position {self.index}")
            name = self._parse_label()
            branch_length = self._parse_branch_length()
            return TreeNode(name=name or None, branch_length=branch_length, children=children)

        name = self._parse_label()
        branch_length = self._parse_branch_length()
        return TreeNode(name=name or None, branch_length=branch_length)

    def _parse_label(self) -> str:
        self._skip_whitespace()
        start = self.index
        while self.index < self.length and self.text[self.index] not in ",:();":
            self.index += 1
        return self.text[start:self.index].strip()

    def _parse_branch_length(self) -> float | None:
        self._skip_whitespace()
        if self._peek() != ":":
            return None
        self.index += 1
        self._skip_whitespace()
        start = self.index
        while self.index < self.length and self.text[self.index] not in ",();":
            self.index += 1
        raw_value = self.text[start:self.index].strip()
        if not raw_value:
            raise ValueError(f"missing branch length in Newick string at position {start}")
        return float(raw_value)

    def _skip_whitespace(self) -> None:
        while self.index < self.length and self.text[self.index].isspace():
            self.index += 1

    def _peek(self) -> str | None:
        if self.index >= self.length:
            return None
        return self.text[self.index]


def loads_newick(text: str) -> PhyloTree:
    """Parse a Newick string into a minimal tree object."""
    return _NewickParser(text).parse()


def load_newick(path: Path) -> PhyloTree:
    """Load a Newick tree from disk."""
    if not path.exists():
        raise FileNotFoundError(f"tree file not found: {path}")
    return loads_newick(path.read_text(encoding="utf-8"))

