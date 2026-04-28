from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.diagnostics.validation import _load_tree
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


def _informative_clade_nodes(tree: PhyloTree, shared_taxa: set[str]) -> dict[frozenset[str], TreeNode]:
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
    try:
        return float(value)
    except ValueError:
        return None


def _load_biophylo_tree(path: Path) -> Tree:
    return Phylo.read(path, detect_tree_format(path))


def _informative_biophylo_clades(tree: Tree, shared_taxa: set[str]) -> dict[frozenset[str], Clade]:
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


def compare_tree_paths(left_path: Path, right_path: Path) -> TreeComparisonReport:
    """Compare two trees over their shared taxa."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
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
    return TreeComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=sorted(shared_taxa),
        left_only_taxa=sorted(left_taxa - right_taxa),
        right_only_taxa=sorted(right_taxa - left_taxa),
        left_informative_clades=len(left_clades),
        right_informative_clades=len(right_clades),
        robinson_foulds_distance=len(symmetric_difference),
        normalized_robinson_foulds=normalized,
    )


def compare_support_values(left_path: Path, right_path: Path) -> SupportComparisonReport:
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
    shared_clade_ids = sorted(left_clades.keys() & right_clades.keys(), key=lambda item: sorted(item))
    return SupportComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=sorted(shared_taxa),
        shared_clades=[
            CladeSupportPair(
                split_id="|".join(sorted(clade_id)),
                left_support=_parse_support(left_clades[clade_id].confidence or left_clades[clade_id].name),
                right_support=_parse_support(right_clades[clade_id].confidence or right_clades[clade_id].name),
            )
            for clade_id in shared_clade_ids
        ],
    )
