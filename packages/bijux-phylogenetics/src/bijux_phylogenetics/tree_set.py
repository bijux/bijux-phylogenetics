from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree as BioTree

from bijux_phylogenetics.compare.topology import _canonical_bipartition, _informative_clades, _unrooted_splits
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import InvalidAlignmentError
from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.io.trees import detect_tree_format


@dataclass(frozen=True, slots=True)
class TreeSetRecord:
    index: int
    tip_count: int
    taxa: list[str]
    rooted_topology_id: str
    unrooted_topology_id: str


@dataclass(slots=True)
class TreeSetReport:
    path: Path
    source_format: str
    tree_count: int
    shared_taxa: list[str]
    taxa_union: list[str]
    rooted_topology_count: int
    unrooted_topology_count: int
    records: list[TreeSetRecord]


@dataclass(frozen=True, slots=True)
class CladeFrequency:
    clade: str
    tree_count: int
    frequency: float


@dataclass(slots=True)
class CladeFrequencyReport:
    path: Path
    tree_count: int
    shared_taxa: list[str]
    clade_frequencies: list[CladeFrequency]


@dataclass(slots=True)
class ConsensusTreeReport:
    path: Path
    tree_count: int
    shared_taxa: list[str]
    consensus_newick: str


@dataclass(frozen=True, slots=True)
class TreeDistancePair:
    left_index: int
    right_index: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float


@dataclass(slots=True)
class TreeDistanceMatrixReport:
    path: Path
    tree_count: int
    shared_taxa: list[str]
    pairs: list[TreeDistancePair]


@dataclass(frozen=True, slots=True)
class TreeTopologyCluster:
    rooted_topology_id: str
    tree_indices: list[int]
    tree_count: int
    representative_index: int


@dataclass(slots=True)
class TreeTopologyClusterReport:
    path: Path
    tree_count: int
    rooted_topology_count: int
    clusters: list[TreeTopologyCluster]


def _shared_taxa(trees: list[PhyloTree]) -> set[str]:
    shared = set(trees[0].tip_names)
    for tree in trees[1:]:
        shared &= set(tree.tip_names)
    return shared


def _taxa_union(trees: list[PhyloTree]) -> set[str]:
    taxa: set[str] = set()
    for tree in trees:
        taxa.update(tree.tip_names)
    return taxa


def _require_tree_set(path: Path) -> tuple[str, list[BioTree], list[PhyloTree]]:
    if not path.exists():
        raise FileNotFoundError(f"tree-set file not found: {path}")
    source_format = detect_tree_format(path)
    bio_trees = list(Phylo.parse(path, source_format))
    if not bio_trees:
        raise InvalidAlignmentError(f"tree set contains no trees: {path}")
    trees = [tree_from_biophylo(tree, source_format=source_format) for tree in bio_trees]
    return source_format, bio_trees, trees


def _format_clade(clade: frozenset[str]) -> str:
    return "|".join(sorted(clade))


def _rooted_topology_id(tree: PhyloTree, shared_taxa: set[str]) -> str:
    return "||".join(sorted(_format_clade(clade) for clade in _informative_clades(tree, shared_taxa)))


def _unrooted_topology_id(tree: PhyloTree, shared_taxa: set[str]) -> str:
    return "||".join(sorted(_format_clade(clade) for clade in _unrooted_splits(tree, shared_taxa)))


def _validate_same_taxa(trees: list[PhyloTree]) -> list[str]:
    first = sorted(trees[0].tip_names)
    for tree in trees[1:]:
        if sorted(tree.tip_names) != first:
            raise InvalidAlignmentError("tree-set analysis requires all trees to share the exact same taxon set")
    return first


def _clade_branch_lengths(tree: PhyloTree, shared_taxa: set[str]) -> dict[frozenset[str], float | None]:
    lengths: dict[frozenset[str], float | None] = {}

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name is not None and node.name in shared_taxa else set()
        taxa: set[str] = set()
        for child in node.children:
            taxa.update(visit(child))
        if 1 < len(taxa) < len(shared_taxa):
            lengths[frozenset(taxa)] = node.branch_length
        return taxa

    visit(tree.root)
    return lengths


def _terminal_branch_lengths(tree: PhyloTree) -> dict[str, float | None]:
    return {name: length for name, length in tree.terminal_branch_lengths() if name is not None}


def _maximal_nested_clades(parent: frozenset[str], clades: set[frozenset[str]]) -> list[frozenset[str]]:
    nested = [clade for clade in clades if clade < parent]
    return sorted(
        [
            clade
            for clade in nested
            if not any(clade < other < parent for other in nested)
        ],
        key=lambda clade: (len(clade), sorted(clade)),
    )


def _mean(values: list[float]) -> float:
    return round(sum(values) / len(values), 15)


def _build_consensus_node(
    taxa: frozenset[str],
    *,
    majority_clades: set[frozenset[str]],
    clade_support: dict[frozenset[str], float],
    clade_lengths: dict[frozenset[str], float],
    terminal_lengths: dict[str, float],
    is_root: bool = False,
) -> TreeNode:
    child_clades = _maximal_nested_clades(taxa, majority_clades)
    covered: set[str] = set()
    children: list[TreeNode] = []
    for child_clade in child_clades:
        covered.update(child_clade)
        children.append(
            _build_consensus_node(
                child_clade,
                majority_clades=majority_clades,
                clade_support=clade_support,
                clade_lengths=clade_lengths,
                terminal_lengths=terminal_lengths,
            )
        )
    for taxon in sorted(taxa - covered):
        children.append(TreeNode(name=taxon, branch_length=terminal_lengths.get(taxon)))
    if len(children) == 1:
        return children[0]
    label = None if is_root else format(clade_support[taxa], ".15g")
    return TreeNode(
        name=label,
        branch_length=None if is_root else clade_lengths.get(taxa),
        children=children,
    )


def _tree_distance(left: PhyloTree, right: PhyloTree, shared_taxa: set[str]) -> tuple[int, float]:
    left_clades = _informative_clades(left, shared_taxa)
    right_clades = _informative_clades(right, shared_taxa)
    symmetric_difference = left_clades.symmetric_difference(right_clades)
    denominator = len(left_clades) + len(right_clades)
    normalized = 0.0 if denominator == 0 else len(symmetric_difference) / denominator
    return len(symmetric_difference), normalized


def load_tree_set(path: Path) -> TreeSetReport:
    """Read a set of trees and summarize their topology diversity over shared taxa."""
    source_format, _, trees = _require_tree_set(path)
    shared_taxa = sorted(_shared_taxa(trees))
    taxa_union = sorted(_taxa_union(trees))
    records = [
        TreeSetRecord(
            index=index,
            tip_count=tree.tip_count,
            taxa=sorted(tree.tip_names),
            rooted_topology_id=_rooted_topology_id(tree, set(shared_taxa)),
            unrooted_topology_id=_unrooted_topology_id(tree, set(shared_taxa)),
        )
        for index, tree in enumerate(trees, start=1)
    ]
    return TreeSetReport(
        path=path,
        source_format=source_format,
        tree_count=len(trees),
        shared_taxa=shared_taxa,
        taxa_union=taxa_union,
        rooted_topology_count=len({record.rooted_topology_id for record in records}),
        unrooted_topology_count=len({record.unrooted_topology_id for record in records}),
        records=records,
    )


def compute_clade_frequency_table(path: Path) -> CladeFrequencyReport:
    """Compute informative clade frequencies across a tree set with a shared taxon set."""
    _, _, trees = _require_tree_set(path)
    shared_taxa = set(_validate_same_taxa(trees))
    counts: dict[str, int] = {}
    for tree in trees:
        for clade in _informative_clades(tree, shared_taxa):
            counts[_format_clade(clade)] = counts.get(_format_clade(clade), 0) + 1
    total = len(trees)
    return CladeFrequencyReport(
        path=path,
        tree_count=total,
        shared_taxa=sorted(shared_taxa),
        clade_frequencies=[
            CladeFrequency(
                clade=clade,
                tree_count=count,
                frequency=round(count / total, 15),
            )
            for clade, count in sorted(counts.items())
        ],
    )


def write_clade_frequency_table(path: Path, report: CladeFrequencyReport) -> Path:
    """Write a clade-frequency table as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["clade", "tree_count", "frequency"], delimiter="\t")
        writer.writeheader()
        for row in report.clade_frequencies:
            writer.writerow(
                {
                    "clade": row.clade,
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                }
            )
    return path


def compute_consensus_tree(path: Path) -> tuple[PhyloTree, ConsensusTreeReport]:
    """Compute a majority-rule consensus tree from a tree set."""
    source_format, _, trees = _require_tree_set(path)
    shared_taxa = _validate_same_taxa(trees)
    universe = frozenset(shared_taxa)
    counts: dict[frozenset[str], int] = {}
    branch_lengths_by_clade: dict[frozenset[str], list[float]] = {}
    terminal_lengths: dict[str, list[float]] = {}
    for tree in trees:
        for clade, length in _clade_branch_lengths(tree, set(shared_taxa)).items():
            counts[clade] = counts.get(clade, 0) + 1
            if length is not None:
                branch_lengths_by_clade.setdefault(clade, []).append(float(length))
        for taxon, length in _terminal_branch_lengths(tree).items():
            if length is not None:
                terminal_lengths.setdefault(taxon, []).append(float(length))
    majority_clades = {
        clade
        for clade, count in counts.items()
        if count / len(trees) > 0.5
    }
    clade_support = {
        clade: round((counts[clade] / len(trees)) * 100.0, 15)
        for clade in majority_clades
    }
    clade_lengths = {
        clade: _mean(lengths)
        for clade, lengths in branch_lengths_by_clade.items()
        if clade in majority_clades and lengths
    }
    terminal_length_means = {
        taxon: _mean(lengths)
        for taxon, lengths in terminal_lengths.items()
        if lengths
    }
    tree = PhyloTree(
        root=_build_consensus_node(
            universe,
            majority_clades=majority_clades,
            clade_support=clade_support,
            clade_lengths=clade_lengths,
            terminal_lengths=terminal_length_means,
            is_root=True,
        ),
        source_format=source_format,
    )
    return tree, ConsensusTreeReport(
        path=path,
        tree_count=len(trees),
        shared_taxa=shared_taxa,
        consensus_newick=dumps_newick(tree),
    )


def write_consensus_tree(path: Path, tree: PhyloTree) -> Path:
    """Write a consensus tree as canonical Newick."""
    return write_newick(path, tree)


def compute_tree_distance_matrix(path: Path) -> TreeDistanceMatrixReport:
    """Compute a pairwise RF-distance matrix across a tree set."""
    _, _, trees = _require_tree_set(path)
    shared_taxa = set(_validate_same_taxa(trees))
    pairs: list[TreeDistancePair] = []
    for left_index, left in enumerate(trees, start=1):
        for right_index, right in enumerate(trees[left_index - 1 :], start=left_index):
            distance, normalized = _tree_distance(left, right, shared_taxa)
            pairs.append(
                TreeDistancePair(
                    left_index=left_index,
                    right_index=right_index,
                    robinson_foulds_distance=distance,
                    normalized_robinson_foulds=normalized,
                )
            )
    return TreeDistanceMatrixReport(
        path=path,
        tree_count=len(trees),
        shared_taxa=sorted(shared_taxa),
        pairs=pairs,
    )


def write_tree_distance_matrix(path: Path, report: TreeDistanceMatrixReport) -> Path:
    """Write a pairwise tree-distance matrix as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "left_index",
                "right_index",
                "robinson_foulds_distance",
                "normalized_robinson_foulds",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.pairs:
            writer.writerow(
                {
                    "left_index": row.left_index,
                    "right_index": row.right_index,
                    "robinson_foulds_distance": row.robinson_foulds_distance,
                    "normalized_robinson_foulds": format(row.normalized_robinson_foulds, ".15g"),
                }
            )
    return path


def cluster_trees_by_topology(path: Path) -> TreeTopologyClusterReport:
    """Cluster trees by identical rooted topology signatures."""
    report = load_tree_set(path)
    clusters_by_id: dict[str, list[int]] = {}
    for record in report.records:
        clusters_by_id.setdefault(record.rooted_topology_id, []).append(record.index)
    clusters = [
        TreeTopologyCluster(
            rooted_topology_id=topology_id,
            tree_indices=indices,
            tree_count=len(indices),
            representative_index=indices[0],
        )
        for topology_id, indices in sorted(
            clusters_by_id.items(),
            key=lambda item: (-len(item[1]), item[1][0]),
        )
    ]
    return TreeTopologyClusterReport(
        path=path,
        tree_count=report.tree_count,
        rooted_topology_count=report.rooted_topology_count,
        clusters=clusters,
    )
