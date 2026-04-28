from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable


@dataclass(slots=True)
class TreeNode:
    """Minimal in-memory tree node."""

    name: str | None = None
    branch_length: float | None = None
    children: list["TreeNode"] = field(default_factory=list)

    def is_leaf(self) -> bool:
        return not self.children

    def iter_nodes(self) -> Iterable["TreeNode"]:
        yield self
        for child in self.children:
            yield from child.iter_nodes()

    def iter_leaves(self) -> Iterable["TreeNode"]:
        for node in self.iter_nodes():
            if node.is_leaf():
                yield node


@dataclass(slots=True)
class PhyloTree:
    """Minimal phylogenetic tree container."""

    root: TreeNode
    source_format: str = "newick"

    def iter_nodes(self) -> Iterable[TreeNode]:
        return self.root.iter_nodes()

    def iter_leaves(self) -> Iterable[TreeNode]:
        return self.root.iter_leaves()

    @property
    def tip_names(self) -> list[str]:
        return [node.name for node in self.iter_leaves() if node.name]

    @property
    def tip_count(self) -> int:
        return sum(1 for _ in self.iter_leaves())

    @property
    def internal_node_count(self) -> int:
        return sum(1 for node in self.iter_nodes() if not node.is_leaf())

    def total_branch_length(self) -> float:
        return sum(node.branch_length or 0.0 for node in self.iter_nodes() if node is not self.root)

    def root_to_tip_lengths(self) -> list[float | None]:
        lengths: list[float | None] = []

        def visit(node: TreeNode, distance: float | None) -> None:
            next_distance = distance
            if node is not self.root:
                if distance is None or node.branch_length is None:
                    next_distance = None
                else:
                    next_distance = distance + node.branch_length
            if node.is_leaf():
                lengths.append(next_distance)
                return
            for child in node.children:
                visit(child, next_distance if node is not self.root else 0.0)

        visit(self.root, 0.0)
        return lengths

