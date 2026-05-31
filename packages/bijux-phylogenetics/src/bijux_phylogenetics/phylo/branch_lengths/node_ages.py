from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    NonUltrametricTreeError,
    UnrootedTreeError,
)


@dataclass(frozen=True, slots=True)
class RootedCladeNodeAgeObservation:
    node_kind: str
    descendant_taxa: list[str]
    age: float


def effectively_rooted(tree: PhyloTree) -> bool:
    return tree.rooted is True or len(tree.root.children) == 2


def rooted_clade_node_age_map(
    tree: PhyloTree,
    *,
    tree_path: Path,
    tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
) -> tuple[float, dict[frozenset[str], RootedCladeNodeAgeObservation]]:
    """Return one ultrametric root age plus one node-age observation per rooted clade."""
    if not effectively_rooted(tree):
        raise UnrootedTreeError(
            "date-aware tree comparison requires rooted trees",
            code="date_aware_tree_comparison_requires_rooted_tree",
            details={"tree_path": str(tree_path), "rooted": tree.rooted},
        )
    tip_depths = _root_to_tip_depths(tree, tree_path=tree_path)
    ultrametric = summarize_ultrametric_tip_depths(tip_depths, tolerance=tolerance)
    if not ultrametric.ultrametric:
        raise NonUltrametricTreeError(
            "date-aware tree comparison requires ultrametric trees",
            code="date_aware_tree_comparison_requires_ultrametric_tree",
            details={
                "tree_path": str(tree_path),
                "minimum_tip_depth": ultrametric.minimum_tip_depth,
                "maximum_tip_depth": ultrametric.maximum_tip_depth,
                "max_tip_depth_deviation": ultrametric.max_tip_depth_deviation,
                "offending_taxa": ultrametric.offending_taxa,
                "tolerance": ultrametric.tolerance,
            },
        )
    root_age = ultrametric.root_age
    clade_ages: dict[frozenset[str], RootedCladeNodeAgeObservation] = {}

    def visit(node: TreeNode, depth: float) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name is not None else set()
        descendant_taxa: set[str] = set()
        for child in node.children:
            if child.branch_length is None:
                raise InvalidBranchLengthError(
                    "date-aware tree comparison requires complete branch lengths",
                    code="date_aware_tree_comparison_requires_complete_branch_lengths",
                    details={"tree_path": str(tree_path), "node_id": child.node_id},
                )
            descendant_taxa.update(visit(child, depth + child.branch_length))
        if len(descendant_taxa) >= 2:
            clade_ages[frozenset(descendant_taxa)] = RootedCladeNodeAgeObservation(
                node_kind="root" if node is tree.root else "internal",
                descendant_taxa=sorted(descendant_taxa),
                age=round(root_age - depth, 15),
            )
        return descendant_taxa

    visit(tree.root, 0.0)
    return root_age, clade_ages


def _root_to_tip_depths(tree: PhyloTree, *, tree_path: Path) -> dict[str, float]:
    tip_depths: dict[str, float] = {}

    def visit(node: TreeNode, depth: float) -> None:
        if node.is_leaf():
            if node.name is not None:
                tip_depths[node.name] = depth
            return
        for child in node.children:
            if child.branch_length is None:
                raise InvalidBranchLengthError(
                    "date-aware tree comparison requires complete branch lengths",
                    code="date_aware_tree_comparison_requires_complete_branch_lengths",
                    details={"tree_path": str(tree_path), "node_id": child.node_id},
                )
            visit(child, depth + child.branch_length)

    visit(tree.root, 0.0)
    return tip_depths
