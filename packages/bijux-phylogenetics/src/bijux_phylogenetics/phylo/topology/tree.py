from __future__ import annotations

from collections.abc import Iterable
from copy import deepcopy
from dataclasses import dataclass, field
import re

_NORMALIZE_TAXON_KEY_PATTERN = re.compile(r"[^0-9A-Za-z._-]+")
_NODE_TOKEN_ESCAPE = str.maketrans(
    {
        "%": "%25",
        "/": "%2F",
        "|": "%7C",
    }
)


def normalize_taxon_key(raw: str) -> str:
    """Normalize a raw taxon label into a stable comparison key."""
    normalized = _NORMALIZE_TAXON_KEY_PATTERN.sub("_", raw.strip())
    return normalized.strip("_")


def _escape_node_token(raw: str) -> str:
    return raw.translate(_NODE_TOKEN_ESCAPE)


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
    """Native in-memory tree node with durable identity and metadata."""

    name: str | None = None
    branch_length: float | None = None
    children: list[TreeNode] = field(default_factory=list)
    node_id: str | None = None
    metadata: dict[str, object] = field(default_factory=dict)
    edge_metadata: dict[str, object] = field(default_factory=dict)
    parent: TreeNode | None = field(default=None, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        self.replace_children(self.children)

    def is_leaf(self) -> bool:
        return not self.children

    def replace_children(self, children: Iterable[TreeNode]) -> None:
        """Replace child links while keeping parent pointers consistent."""
        for child in self.children:
            if child.parent is self:
                child.parent = None
        self.children = list(children)
        for child in self.children:
            child.parent = self

    def append_child(self, child: TreeNode) -> None:
        child.parent = self
        self.children.append(child)

    def extend_children(self, children: Iterable[TreeNode]) -> None:
        for child in children:
            self.append_child(child)

    def pop_child(self, index: int = -1) -> TreeNode:
        child = self.children.pop(index)
        if child.parent is self:
            child.parent = None
        return child

    def iter_nodes(self, *, order: str = "preorder") -> Iterable[TreeNode]:
        if order == "preorder":
            yield self
            for child in self.children:
                yield from child.iter_nodes(order=order)
            return
        if order == "postorder":
            for child in self.children:
                yield from child.iter_nodes(order=order)
            yield self
            return
        raise ValueError("tree traversal order must be 'preorder' or 'postorder'")

    def iter_leaves(self) -> Iterable[TreeNode]:
        for node in self.iter_nodes(order="preorder"):
            if node.is_leaf():
                yield node

    def iter_internal_nodes(self, *, order: str = "preorder") -> Iterable[TreeNode]:
        for node in self.iter_nodes(order=order):
            if not node.is_leaf():
                yield node

    def iter_edges(self) -> Iterable[tuple[TreeNode, TreeNode]]:
        for child in self.children:
            yield self, child
            yield from child.iter_edges()

    @property
    def taxon_label(self) -> TaxonLabel | None:
        if self.name is None:
            return None
        return TaxonLabel.from_raw(self.name)

    @property
    def descendant_taxa(self) -> list[str]:
        return descendant_taxa(self)

    def copy(self) -> TreeNode:
        """Deep-copy a node subtree without retaining shared mutable state."""
        return TreeNode(
            name=self.name,
            branch_length=self.branch_length,
            children=[child.copy() for child in self.children],
            node_id=self.node_id,
            metadata=deepcopy(self.metadata),
            edge_metadata=deepcopy(self.edge_metadata),
        )


def descendant_taxa(node: TreeNode) -> list[str]:
    """Return the sorted descendant-tip labels for one node."""
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(descendant_taxa(child))
    return sorted(taxa)


def stable_node_label(node: TreeNode) -> str:
    """Return one durable, human-readable node label token."""
    if node.is_leaf():
        return f"taxon:{_escape_node_token(node.name or '<unnamed-tip>')}"
    taxa = descendant_taxa(node)
    if taxa:
        clade = "|".join(_escape_node_token(taxon) for taxon in taxa)
    else:
        clade = "<unnamed-clade>"
    if node.name is None:
        return f"clade:{clade}"
    return f"clade:{clade};label:{_escape_node_token(node.name)}"


@dataclass(slots=True)
class PhyloTree:
    """Native phylogenetic tree container with stable node identity."""

    root: TreeNode
    source_format: str = "newick"
    rooted: bool | None = None
    metadata: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.refresh()

    def refresh(self) -> PhyloTree:
        """Rebuild parent pointers and stable node identifiers."""
        self._refresh_node(self.root, parent=None, prefix=None)
        return self

    def _refresh_node(
        self,
        node: TreeNode,
        *,
        parent: TreeNode | None,
        prefix: str | None,
    ) -> None:
        node.parent = parent
        segment = stable_node_label(node)
        node.node_id = f"root:{segment}" if prefix is None else f"{prefix}/{segment}"
        for child in node.children:
            self._refresh_node(child, parent=node, prefix=node.node_id)

    def iter_nodes(self, *, order: str = "preorder") -> Iterable[TreeNode]:
        return self.root.iter_nodes(order=order)

    def iter_leaves(self) -> Iterable[TreeNode]:
        return self.root.iter_leaves()

    def iter_internal_nodes(self, *, order: str = "preorder") -> Iterable[TreeNode]:
        return self.root.iter_internal_nodes(order=order)

    def iter_edges(self) -> Iterable[tuple[TreeNode, TreeNode]]:
        return self.root.iter_edges()

    def node_by_id(self, node_id: str) -> TreeNode:
        for node in self.iter_nodes(order="preorder"):
            if node.node_id == node_id:
                return node
        raise KeyError(f"tree does not contain node_id '{node_id}'")

    def copy(self) -> PhyloTree:
        """Deep-copy a tree without shared node, edge, or metadata state."""
        return PhyloTree(
            root=self.root.copy(),
            source_format=self.source_format,
            rooted=self.rooted,
            metadata=deepcopy(self.metadata),
        )

    def validation_errors(self) -> list[str]:
        """Collect structural validation errors for the in-memory tree model."""
        errors: list[str] = []
        seen: set[int] = set()
        stack: set[int] = set()
        node_ids: set[str] = set()

        def visit(node: TreeNode, *, expected_parent: TreeNode | None) -> None:
            identifier = id(node)
            if node.parent is not expected_parent:
                errors.append(
                    f"node '{node.node_id or '<unassigned>'}' has an inconsistent parent pointer"
                )
            if node.node_id is None:
                errors.append("tree contains a node without a stable node_id")
            elif node.node_id in node_ids:
                errors.append(f"tree contains duplicate node_id '{node.node_id}'")
            else:
                node_ids.add(node.node_id)
            if identifier in stack:
                errors.append(
                    f"tree contains a cycle at node '{node.node_id or '<unassigned>'}'"
                )
                return
            if identifier in seen:
                errors.append(
                    f"tree contains duplicate parentage for node '{node.node_id or '<unassigned>'}'"
                )
                return
            seen.add(identifier)
            stack.add(identifier)
            for child in node.children:
                if child is node:
                    errors.append(
                        f"node '{node.node_id or '<unassigned>'}' references itself as a child"
                    )
                    continue
                visit(child, expected_parent=node)
            stack.remove(identifier)

        if self.root.parent is not None:
            errors.append("tree root must not have a parent")
        visit(self.root, expected_parent=None)
        if not any(True for _ in self.iter_leaves()):
            errors.append("tree must contain at least one terminal node")
        return errors

    def validate(self) -> None:
        """Raise if the in-memory tree structure is internally inconsistent."""
        errors = self.validation_errors()
        if errors:
            raise ValueError("; ".join(errors))

    def to_newick(self) -> str:
        """Serialize the tree through the repository's canonical Newick writer."""
        from bijux_phylogenetics.io.newick import dumps_newick

        return dumps_newick(self)

    @classmethod
    def from_newick(cls, text: str) -> PhyloTree:
        """Parse one Newick string through the repository's current tree loader."""
        from bijux_phylogenetics.io.newick import loads_newick

        return loads_newick(text)

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
        return sum(1 for _ in self.iter_internal_nodes())

    def total_branch_length(self) -> float:
        return sum(
            node.branch_length or 0.0
            for node in self.iter_nodes(order="preorder")
            if node is not self.root
        )

    def branch_lengths(self) -> list[float | None]:
        return [
            node.branch_length
            for node in self.iter_nodes(order="preorder")
            if node is not self.root
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
