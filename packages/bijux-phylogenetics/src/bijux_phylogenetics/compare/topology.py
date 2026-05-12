from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from bijux_phylogenetics.core.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.diagnostics.validation import _load_tree
from bijux_phylogenetics.io.iqtree_support import parse_iqtree_branch_support_label
from bijux_phylogenetics.io.trees import detect_tree_format


@dataclass(slots=True)
class TreeComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    left_informative_clades: int
    right_informative_clades: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    topology_equal: bool
    same_unrooted_topology: bool
    same_taxa_different_rooting: bool
    same_topology_different_branch_lengths: bool


@dataclass(slots=True)
class SharedTaxaPruningReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]


@dataclass(slots=True)
class CladeSupportPair:
    split_id: str
    left_support: float | None
    right_support: float | None


@dataclass(slots=True)
class SupportComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    shared_clades: list[CladeSupportPair]


@dataclass(slots=True)
class BranchLengthPair:
    split_id: str
    left_length: float | None
    right_length: float | None
    delta: float | None
    ratio: float | None


@dataclass(slots=True)
class BranchLengthComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    shared_splits: list[BranchLengthPair]


@dataclass(slots=True)
class CladeSetComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    shared_clades: list[str]
    left_only_clades: list[str]
    right_only_clades: list[str]


@dataclass(slots=True)
class CladeChangeReport:
    left_path: Path
    right_path: Path
    lost_clades: list[str]
    gained_clades: list[str]


@dataclass(slots=True)
class InMemoryTopologyComparison:
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    left_informative_clades: int
    right_informative_clades: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    topology_equal: bool
    same_unrooted_topology: bool
    same_taxa_different_rooting: bool


@dataclass(slots=True)
class InMemoryBranchLengthComparison:
    shared_taxa: list[str]
    shared_splits: list[BranchLengthPair]


def _informative_clades(tree: PhyloTree, shared_taxa: set[str]) -> set[frozenset[str]]:
    clades: set[frozenset[str]] = set()

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name in shared_taxa else set()

        taxa: set[str] = set()
        for child in node.children:
            taxa.update(visit(child))

        if 1 < len(taxa) < len(shared_taxa):
            clades.add(frozenset(taxa))
        return taxa

    visit(tree.root)
    return clades


def _format_clade_set(clades: set[frozenset[str]]) -> list[str]:
    return sorted("|".join(sorted(clade)) for clade in clades)


def _canonical_bipartition(taxa: set[str], universe: set[str]) -> frozenset[str]:
    complement = universe - taxa
    left = sorted(taxa)
    right = sorted(complement)
    if (len(left), left) <= (len(right), right):
        return frozenset(taxa)
    return frozenset(complement)


def _unrooted_splits(tree: PhyloTree, shared_taxa: set[str]) -> set[frozenset[str]]:
    splits: set[frozenset[str]] = set()

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name in shared_taxa else set()

        taxa: set[str] = set()
        for child in node.children:
            taxa.update(visit(child))
        if node is not tree.root and 1 < len(taxa) < len(shared_taxa) - 1:
            splits.add(_canonical_bipartition(taxa, shared_taxa))
        return taxa

    visit(tree.root)
    return splits


def _informative_clade_nodes(
    tree: PhyloTree, shared_taxa: set[str]
) -> dict[frozenset[str], TreeNode]:
    clades: dict[frozenset[str], TreeNode] = {}

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name in shared_taxa else set()

        taxa: set[str] = set()
        for child in node.children:
            taxa.update(visit(child))

        if 1 < len(taxa) < len(shared_taxa):
            clades[frozenset(taxa)] = node
        return taxa

    visit(tree.root)
    return clades


def _parse_support(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        parsed = parse_iqtree_branch_support_label(value)
        if parsed is not None:
            return (
                parsed.ufboot_support
                if parsed.ufboot_support is not None
                else parsed.sh_alrt_support
            )
    try:
        return float(value)
    except ValueError:
        return None


def _load_biophylo_tree(path: Path) -> Tree:
    return Phylo.read(path, detect_tree_format(path))


def _compare_tree_objects(
    left: PhyloTree, right: PhyloTree
) -> InMemoryTopologyComparison:
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa = left_taxa & right_taxa
    if len(shared_taxa) < 2:
        raise ValueError("tree comparison requires at least two shared taxa")

    left_clades = _informative_clades(left, shared_taxa)
    right_clades = _informative_clades(right, shared_taxa)
    symmetric_difference = left_clades.symmetric_difference(right_clades)
    denominator = len(left_clades) + len(right_clades)
    normalized = 0.0 if denominator == 0 else len(symmetric_difference) / denominator
    topology_equal = len(symmetric_difference) == 0
    same_unrooted_topology = _unrooted_splits(left, shared_taxa) == _unrooted_splits(
        right, shared_taxa
    )
    same_taxa_different_rooting = (
        left_taxa == right_taxa and same_unrooted_topology and not topology_equal
    )
    return InMemoryTopologyComparison(
        shared_taxa=sorted(shared_taxa),
        left_only_taxa=sorted(left_taxa - right_taxa),
        right_only_taxa=sorted(right_taxa - left_taxa),
        left_informative_clades=len(left_clades),
        right_informative_clades=len(right_clades),
        robinson_foulds_distance=len(symmetric_difference),
        normalized_robinson_foulds=normalized,
        topology_equal=topology_equal,
        same_unrooted_topology=same_unrooted_topology,
        same_taxa_different_rooting=same_taxa_different_rooting,
    )


def _compare_branch_lengths_for_trees(
    left: PhyloTree, right: PhyloTree
) -> InMemoryBranchLengthComparison:
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa = left_taxa & right_taxa
    if len(shared_taxa) < 2:
        raise ValueError("branch-length comparison requires at least two shared taxa")

    left_clades = _informative_clade_nodes(left, shared_taxa)
    right_clades = _informative_clade_nodes(right, shared_taxa)
    shared_clade_ids = sorted(
        left_clades.keys() & right_clades.keys(), key=lambda item: sorted(item)
    )
    pairs: list[BranchLengthPair] = []
    for clade_id in shared_clade_ids:
        left_length = left_clades[clade_id].branch_length
        right_length = right_clades[clade_id].branch_length
        delta = (
            None
            if left_length is None or right_length is None
            else right_length - left_length
        )
        ratio = None
        if left_length not in {None, 0} and right_length is not None:
            ratio = right_length / left_length
        pairs.append(
            BranchLengthPair(
                split_id="|".join(sorted(clade_id)),
                left_length=left_length,
                right_length=right_length,
                delta=delta,
                ratio=ratio,
            )
        )
    return InMemoryBranchLengthComparison(
        shared_taxa=sorted(shared_taxa), shared_splits=pairs
    )


def _informative_biophylo_clades(
    tree: Tree, shared_taxa: set[str]
) -> dict[frozenset[str], Clade]:
    clades: dict[frozenset[str], Clade] = {}

    def visit(clade: Clade) -> set[str]:
        if not clade.clades:
            return {clade.name} if clade.name in shared_taxa else set()

        taxa: set[str] = set()
        for child in clade.clades:
            taxa.update(visit(child))

        if 1 < len(taxa) < len(shared_taxa):
            clades[frozenset(taxa)] = clade
        return taxa

    visit(tree.root)
    return clades


def prune_trees_to_shared_taxa(
    left_path: Path,
    right_path: Path,
) -> tuple[PhyloTree, PhyloTree, SharedTaxaPruningReport]:
    """Prune two trees to the exact shared taxon set."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa = sorted(left_taxa & right_taxa)
    if len(shared_taxa) < 2:
        raise ValueError("shared-taxon pruning requires at least two shared taxa")

    pruned_left, _ = prune_tree_to_requested_taxa(left_path, shared_taxa)
    pruned_right, _ = prune_tree_to_requested_taxa(right_path, shared_taxa)
    return (
        pruned_left,
        pruned_right,
        SharedTaxaPruningReport(
            left_path=left_path,
            right_path=right_path,
            shared_taxa=shared_taxa,
            left_only_taxa=sorted(left_taxa - right_taxa),
            right_only_taxa=sorted(right_taxa - left_taxa),
        ),
    )


def compare_clade_sets(left_path: Path, right_path: Path) -> CladeSetComparisonReport:
    """Compare rooted informative clade sets across two trees."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    left_taxa = set(left.tip_names)
    right_taxa = set(right.tip_names)
    shared_taxa = left_taxa & right_taxa
    if len(shared_taxa) < 2:
        raise ValueError("clade comparison requires at least two shared taxa")

    left_clades = _informative_clades(left, shared_taxa)
    right_clades = _informative_clades(right, shared_taxa)
    return CladeSetComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=sorted(shared_taxa),
        shared_clades=_format_clade_set(left_clades & right_clades),
        left_only_clades=_format_clade_set(left_clades - right_clades),
        right_only_clades=_format_clade_set(right_clades - left_clades),
    )


def detect_clade_changes(left_path: Path, right_path: Path) -> CladeChangeReport:
    """Report clades lost from the left tree and gained in the right tree."""
    report = compare_clade_sets(left_path, right_path)
    return CladeChangeReport(
        left_path=left_path,
        right_path=right_path,
        lost_clades=report.left_only_clades,
        gained_clades=report.right_only_clades,
    )


def compare_tree_paths(left_path: Path, right_path: Path) -> TreeComparisonReport:
    """Compare two trees over their shared taxa."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    comparison = _compare_tree_objects(left, right)
    branch_report = _compare_branch_lengths_for_trees(left, right)
    same_topology_different_branch_lengths = comparison.topology_equal and any(
        row.left_length != row.right_length for row in branch_report.shared_splits
    )
    return TreeComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=comparison.shared_taxa,
        left_only_taxa=comparison.left_only_taxa,
        right_only_taxa=comparison.right_only_taxa,
        left_informative_clades=comparison.left_informative_clades,
        right_informative_clades=comparison.right_informative_clades,
        robinson_foulds_distance=comparison.robinson_foulds_distance,
        normalized_robinson_foulds=comparison.normalized_robinson_foulds,
        topology_equal=comparison.topology_equal,
        same_unrooted_topology=comparison.same_unrooted_topology,
        same_taxa_different_rooting=comparison.same_taxa_different_rooting,
        same_topology_different_branch_lengths=same_topology_different_branch_lengths,
    )


def compare_support_values(
    left_path: Path, right_path: Path
) -> SupportComparisonReport:
    """Compare internal clade support values across two trees with shared taxa."""
    left = _load_biophylo_tree(left_path)
    right = _load_biophylo_tree(right_path)
    left_taxa = {clade.name for clade in left.get_terminals() if clade.name}
    right_taxa = {clade.name for clade in right.get_terminals() if clade.name}
    shared_taxa = left_taxa & right_taxa
    if len(shared_taxa) < 2:
        raise ValueError("support comparison requires at least two shared taxa")

    left_clades = _informative_biophylo_clades(left, shared_taxa)
    right_clades = _informative_biophylo_clades(right, shared_taxa)
    shared_clade_ids = sorted(
        left_clades.keys() & right_clades.keys(), key=lambda item: sorted(item)
    )
    return SupportComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=sorted(shared_taxa),
        shared_clades=[
            CladeSupportPair(
                split_id="|".join(sorted(clade_id)),
                left_support=_parse_support(
                    left_clades[clade_id].confidence or left_clades[clade_id].name
                ),
                right_support=_parse_support(
                    right_clades[clade_id].confidence or right_clades[clade_id].name
                ),
            )
            for clade_id in shared_clade_ids
        ],
    )


def compare_branch_lengths(
    left_path: Path, right_path: Path
) -> BranchLengthComparisonReport:
    """Compare branch lengths across matching informative splits."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    comparison = _compare_branch_lengths_for_trees(left, right)
    return BranchLengthComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=comparison.shared_taxa,
        shared_splits=comparison.shared_splits,
    )


def write_tree_comparison_table(path: Path, left_path: Path, right_path: Path) -> Path:
    """Write a flat TSV table covering the compared clade and split surfaces."""
    clades = compare_clade_sets(left_path, right_path)
    support = compare_support_values(left_path, right_path)
    branch_lengths = compare_branch_lengths(left_path, right_path)
    support_by_id = {row.split_id: row for row in support.shared_clades}
    branch_by_id = {row.split_id: row for row in branch_lengths.shared_splits}
    all_split_ids = sorted(
        set(clades.shared_clades)
        | set(clades.left_only_clades)
        | set(clades.right_only_clades)
        | set(support_by_id)
        | set(branch_by_id)
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "split_id",
                "comparison_status",
                "shared_clade",
                "left_support",
                "right_support",
                "left_length",
                "right_length",
                "length_delta",
                "length_ratio",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for split_id in all_split_ids:
            support_row = support_by_id.get(split_id)
            branch_row = branch_by_id.get(split_id)
            if split_id in clades.shared_clades:
                status = "shared"
            elif split_id in clades.left_only_clades:
                status = "left_only"
            else:
                status = "right_only"
            writer.writerow(
                {
                    "split_id": split_id,
                    "comparison_status": status,
                    "shared_clade": str(split_id in clades.shared_clades).lower(),
                    "left_support": ""
                    if support_row is None or support_row.left_support is None
                    else support_row.left_support,
                    "right_support": ""
                    if support_row is None or support_row.right_support is None
                    else support_row.right_support,
                    "left_length": ""
                    if branch_row is None or branch_row.left_length is None
                    else branch_row.left_length,
                    "right_length": ""
                    if branch_row is None or branch_row.right_length is None
                    else branch_row.right_length,
                    "length_delta": ""
                    if branch_row is None or branch_row.delta is None
                    else branch_row.delta,
                    "length_ratio": ""
                    if branch_row is None or branch_row.ratio is None
                    else branch_row.ratio,
                }
            )
    return path
