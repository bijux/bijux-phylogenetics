from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import _load_tree
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.clades import (
    canonical_clade_id,
    split_sort_key,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    NonUltrametricTreeError,
    UnrootedTreeError,
)

from .comparison import _build_tree_comparison_report
from .models import CladeAgeComparisonRow, DateAwareTreeComparisonReport


@dataclass(frozen=True, slots=True)
class _CladeAgeObservation:
    node_kind: str
    descendant_taxa: list[str]
    age: float


def _effectively_rooted(tree: PhyloTree) -> bool:
    return tree.rooted is True or len(tree.root.children) == 2


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


def _clade_age_map(
    tree: PhyloTree,
    *,
    tree_path: Path,
    tolerance: float,
) -> tuple[float, dict[frozenset[str], _CladeAgeObservation]]:
    if not _effectively_rooted(tree):
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
    clade_ages: dict[frozenset[str], _CladeAgeObservation] = {}

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
            clade_ages[frozenset(descendant_taxa)] = _CladeAgeObservation(
                node_kind="root" if node is tree.root else "internal",
                descendant_taxa=sorted(descendant_taxa),
                age=round(root_age - depth, 15),
            )
        return descendant_taxa

    visit(tree.root, 0.0)
    return root_age, clade_ages


def _comparison_scope_trees(
    left_path: Path,
    right_path: Path,
    *,
    comparison_shared_taxa: list[str],
    left_only_taxa: list[str],
    right_only_taxa: list[str],
    taxon_overlap_policy: str,
) -> tuple[PhyloTree, PhyloTree, str]:
    if taxon_overlap_policy == "prune-to-shared" and (left_only_taxa or right_only_taxa):
        pruned_left, _ = prune_tree_to_requested_taxa(left_path, comparison_shared_taxa)
        pruned_right, _ = prune_tree_to_requested_taxa(
            right_path, comparison_shared_taxa
        )
        return pruned_left, pruned_right, "pruned-to-shared-taxa"
    return _load_tree(left_path), _load_tree(right_path), "full-taxa"


def compare_clade_ages(
    left_path: Path,
    right_path: Path,
    *,
    taxon_overlap_policy: str = "prune-to-shared",
    tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
) -> DateAwareTreeComparisonReport:
    """Compare matched rooted clades by node age on one shared taxon scope."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    topology = _build_tree_comparison_report(
        left_path,
        right_path,
        left,
        right,
        rf_mode="rooted",
        taxon_overlap_policy=taxon_overlap_policy,
    )
    comparison_left, comparison_right, comparison_scope = _comparison_scope_trees(
        left_path,
        right_path,
        comparison_shared_taxa=topology.shared_taxa,
        left_only_taxa=topology.left_only_taxa,
        right_only_taxa=topology.right_only_taxa,
        taxon_overlap_policy=taxon_overlap_policy,
    )
    left_root_age, left_clades = _clade_age_map(
        comparison_left,
        tree_path=left_path,
        tolerance=tolerance,
    )
    right_root_age, right_clades = _clade_age_map(
        comparison_right,
        tree_path=right_path,
        tolerance=tolerance,
    )
    matched_signatures = sorted(left_clades.keys() & right_clades.keys(), key=split_sort_key)
    age_differences = [
        round(right_clades[signature].age - left_clades[signature].age, 15)
        for signature in matched_signatures
    ]
    age_rmse = 0.0
    if age_differences:
        age_rmse = round(
            math.sqrt(
                sum(difference * difference for difference in age_differences)
                / len(age_differences)
            ),
            15,
        )
    unstable_age_threshold = age_rmse
    clade_rows = [
        CladeAgeComparisonRow(
            clade_id=canonical_clade_id(signature),
            node_kind=left_clades[signature].node_kind,
            taxon_count=len(signature),
            descendant_taxa=left_clades[signature].descendant_taxa,
            left_age=left_clades[signature].age,
            right_age=right_clades[signature].age,
            age_difference=round(
                right_clades[signature].age - left_clades[signature].age, 15
            ),
            absolute_age_difference=round(
                abs(right_clades[signature].age - left_clades[signature].age), 15
            ),
            unstable_age=(
                abs(right_clades[signature].age - left_clades[signature].age)
                > unstable_age_threshold
            ),
        )
        for signature in matched_signatures
    ]
    absolute_differences = [row.absolute_age_difference for row in clade_rows]
    return DateAwareTreeComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=topology.shared_taxa,
        left_only_taxa=topology.left_only_taxa,
        right_only_taxa=topology.right_only_taxa,
        taxon_overlap_policy=topology.taxon_overlap_policy,
        comparison_scope=comparison_scope,
        left_root_age=left_root_age,
        right_root_age=right_root_age,
        age_rmse=age_rmse,
        mean_absolute_age_difference=(
            0.0
            if not absolute_differences
            else round(sum(absolute_differences) / len(absolute_differences), 15)
        ),
        max_absolute_age_difference=max(absolute_differences, default=0.0),
        unstable_age_threshold=unstable_age_threshold,
        unstable_clade_count=sum(1 for row in clade_rows if row.unstable_age),
        matched_clade_count=len(clade_rows),
        topology=topology,
        clade_rows=clade_rows,
    )
