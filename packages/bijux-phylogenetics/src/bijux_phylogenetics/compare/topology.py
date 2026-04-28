from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.diagnostics.validation import _load_tree


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

