from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
import re

_NORMALIZE_TAXON_KEY_PATTERN = re.compile(r"[^0-9A-Za-z._-]+")


def normalize_taxon_key(raw: str) -> str:
    """Normalize a raw taxon label into a stable comparison key."""
    normalized = _NORMALIZE_TAXON_KEY_PATTERN.sub("_", raw.strip())
    return normalized.strip("_")


@dataclass(frozen=True, slots=True)
class TaxonLabel:
    """Taxon identity preserving the raw label and its normalized key."""

    raw: str
    key: str

    @classmethod
    def from_raw(cls, raw: str) -> TaxonLabel:
        return cls(raw=raw, key=normalize_taxon_key(raw))


@dataclass(slots=True)
class TreeNode:
    """Minimal in-memory tree node."""

    name: str | None = None
    branch_length: float | None = None
    children: list[TreeNode] = field(default_factory=list)

    def is_leaf(self) -> bool:
        return not self.children

    def iter_nodes(self) -> Iterable[TreeNode]:
        yield self
        for child in self.children:
            yield from child.iter_nodes()

    def iter_leaves(self) -> Iterable[TreeNode]:
        for node in self.iter_nodes():
            if node.is_leaf():
                yield node

    @property
    def taxon_label(self) -> TaxonLabel | None:
        if self.name is None:
            return None
        return TaxonLabel.from_raw(self.name)


@dataclass(slots=True)
class PhyloTree:
    """Minimal phylogenetic tree container."""

    root: TreeNode
    source_format: str = "newick"
    rooted: bool | None = None

    def iter_nodes(self) -> Iterable[TreeNode]:
        return self.root.iter_nodes()

    def iter_leaves(self) -> Iterable[TreeNode]:
        return self.root.iter_leaves()

    @property
    def tip_names(self) -> list[str]:
        return [node.name for node in self.iter_leaves() if node.name]

    @property
    def tip_taxa(self) -> list[TaxonLabel]:
        return [
            taxon
            for taxon in (node.taxon_label for node in self.iter_leaves())
            if taxon is not None
        ]

    @property
    def tip_count(self) -> int:
        return sum(1 for _ in self.iter_leaves())

    @property
    def internal_node_count(self) -> int:
        return sum(1 for node in self.iter_nodes() if not node.is_leaf())

    def total_branch_length(self) -> float:
        return sum(
            node.branch_length or 0.0
            for node in self.iter_nodes()
            if node is not self.root
        )

    def branch_lengths(self) -> list[float | None]:
        return [
            node.branch_length for node in self.iter_nodes() if node is not self.root
        ]

    def terminal_branch_lengths(self) -> list[tuple[str, float | None]]:
        return [
            (node.name, node.branch_length)
            for node in self.iter_leaves()
            if node.name is not None
        ]

    def root_to_tip_pairs(self) -> list[tuple[str | None, float | None]]:
        pairs: list[tuple[str | None, float | None]] = []

        def visit(node: TreeNode, distance: float | None) -> None:
            next_distance = distance
            if node is not self.root:
                if distance is None or node.branch_length is None:
                    next_distance = None
                else:
                    next_distance = distance + node.branch_length
            if node.is_leaf():
                pairs.append((node.name, next_distance))
                return
            for child in node.children:
                visit(child, next_distance if node is not self.root else 0.0)

        visit(self.root, 0.0)
        return pairs

    def root_to_tip_lengths(self) -> list[float | None]:
        return [distance for _, distance in self.root_to_tip_pairs()]
