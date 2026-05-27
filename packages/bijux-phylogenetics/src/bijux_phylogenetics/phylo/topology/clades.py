from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json

from bijux_phylogenetics.io.iqtree_support import parse_iqtree_branch_support_label
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

RobinsonFouldsMode = str


@dataclass(frozen=True, slots=True)
class RobinsonFouldsMetrics:
    """Canonical split-set comparison summary for one RF-style tree comparison."""

    rf_mode: str
    left_signatures: frozenset[frozenset[str]]
    right_signatures: frozenset[frozenset[str]]
    shared_signatures: frozenset[frozenset[str]]
    left_only_signatures: frozenset[frozenset[str]]
    right_only_signatures: frozenset[frozenset[str]]
    distance: int
    normalized_distance: float

    @property
    def left_count(self) -> int:
        return len(self.left_signatures)

    @property
    def right_count(self) -> int:
        return len(self.right_signatures)


def canonical_clade_id(signature: frozenset[str]) -> str:
    """Render one descendant-taxon signature into a durable clade identifier."""
    return "|".join(sorted(signature))


def canonical_bipartition(
    descendant_taxa: set[str], universe: set[str]
) -> frozenset[str]:
    """Normalize one unrooted split so child order and side choice do not matter."""
    complement = universe - descendant_taxa
    left = sorted(descendant_taxa)
    right = sorted(complement)
    if (len(left), left) <= (len(right), right):
        return frozenset(descendant_taxa)
    return frozenset(complement)


def informative_rooted_clades(
    tree: PhyloTree,
    shared_taxa: set[str] | None = None,
    *,
    include_root: bool = False,
) -> set[frozenset[str]]:
    """Extract rooted non-singleton clades over one exact or shared taxon scope."""
    clades = informative_rooted_clade_nodes(
        tree,
        shared_taxa,
        include_root=include_root,
    )
    return set(clades)


def informative_rooted_clade_nodes(
    tree: PhyloTree,
    shared_taxa: set[str] | None = None,
    *,
    include_root: bool = False,
) -> dict[frozenset[str], TreeNode]:
    """Map rooted descendant-taxon signatures to the native node that realizes them."""
    taxon_scope = set(tree.tip_names) if shared_taxa is None else set(shared_taxa)
    clades: dict[frozenset[str], TreeNode] = {}
    if len(taxon_scope) < 2:
        return clades

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name in taxon_scope else set()

        descendant_taxa: set[str] = set()
        for child in node.children:
            descendant_taxa.update(visit(child))
        if not descendant_taxa:
            return set()
        if node is tree.root:
            if include_root and len(descendant_taxa) == len(taxon_scope):
                clades[frozenset(descendant_taxa)] = node
            return descendant_taxa
        if 1 < len(descendant_taxa) < len(taxon_scope):
            clades[frozenset(descendant_taxa)] = node
        return descendant_taxa

    visit(tree.root)
    return clades


def informative_unrooted_splits(
    tree: PhyloTree,
    shared_taxa: set[str] | None = None,
) -> set[frozenset[str]]:
    """Extract canonical unrooted bipartitions over one exact or shared taxon scope."""
    taxon_scope = set(tree.tip_names) if shared_taxa is None else set(shared_taxa)
    if len(taxon_scope) < 4:
        return set()
    splits: set[frozenset[str]] = set()

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name in taxon_scope else set()

        descendant_taxa: set[str] = set()
        for child in node.children:
            descendant_taxa.update(visit(child))
        if node is not tree.root and 1 < len(descendant_taxa) < len(taxon_scope) - 1:
            splits.add(canonical_bipartition(descendant_taxa, taxon_scope))
        return descendant_taxa

    visit(tree.root)
    return splits


def tree_has_polytomy(tree: PhyloTree) -> bool:
    """Report whether any node has more than two child branches."""
    return any(len(node.children) > 2 for node in tree.iter_nodes())


def node_support_value(node: TreeNode) -> float | None:
    """Resolve one native node support label using the same interpretation policy everywhere."""
    confidence = node.metadata.get("confidence")
    if confidence is not None:
        return float(confidence)
    if node.name is None:
        return None
    parsed = parse_iqtree_branch_support_label(node.name)
    if parsed is not None:
        return (
            parsed.ufboot_support
            if parsed.ufboot_support is not None
            else parsed.sh_alrt_support
        )
    try:
        return float(node.name)
    except ValueError:
        return None


def robinson_foulds_metrics(
    left: PhyloTree,
    right: PhyloTree,
    shared_taxa: set[str],
    *,
    rf_mode: RobinsonFouldsMode,
) -> RobinsonFouldsMetrics:
    """Compare rooted clades or unrooted splits in one canonical native core."""
    if rf_mode == "rooted":
        left_signatures = informative_rooted_clades(left, shared_taxa)
        right_signatures = informative_rooted_clades(right, shared_taxa)
    elif rf_mode == "unrooted":
        left_signatures = informative_unrooted_splits(left, shared_taxa)
        right_signatures = informative_unrooted_splits(right, shared_taxa)
    else:
        raise ValueError(
            f"rf_mode must be one of {{'rooted', 'unrooted'}}, got {rf_mode!r}"
        )
    shared_signatures = frozenset(left_signatures & right_signatures)
    left_only_signatures = frozenset(left_signatures - right_signatures)
    right_only_signatures = frozenset(right_signatures - left_signatures)
    distance = len(left_only_signatures) + len(right_only_signatures)
    denominator = len(left_signatures) + len(right_signatures)
    normalized = 0.0 if denominator == 0 else distance / denominator
    return RobinsonFouldsMetrics(
        rf_mode=rf_mode,
        left_signatures=frozenset(left_signatures),
        right_signatures=frozenset(right_signatures),
        shared_signatures=shared_signatures,
        left_only_signatures=left_only_signatures,
        right_only_signatures=right_only_signatures,
        distance=distance,
        normalized_distance=normalized,
    )


def split_sort_key(signature: frozenset[str]) -> tuple[int, tuple[str, ...]]:
    """Sort clades and splits deterministically for ledgers and reports."""
    ordered_taxa = tuple(sorted(signature))
    return (len(ordered_taxa), ordered_taxa)


def rooted_topology_signature_ids(
    tree: PhyloTree,
    shared_taxa: set[str] | None = None,
) -> tuple[str, ...]:
    """Render one rooted tree topology as sorted informative clade identifiers."""
    clades = informative_rooted_clades(tree, shared_taxa)
    return tuple(
        canonical_clade_id(signature)
        for signature in sorted(clades, key=split_sort_key)
    )


def rooted_topology_fingerprint(
    tree: PhyloTree,
    shared_taxa: set[str] | None = None,
) -> str:
    """Hash one rooted topology independent of branch lengths and child order."""
    taxon_scope = sorted(tree.tip_names) if shared_taxa is None else sorted(shared_taxa)
    payload = {
        "taxa": taxon_scope,
        "informative_rooted_clades": rooted_topology_signature_ids(tree, shared_taxa),
    }
    return hashlib.sha256(
        json.dumps(payload, separators=(",", ":"), sort_keys=True).encode("utf-8")
    ).hexdigest()
