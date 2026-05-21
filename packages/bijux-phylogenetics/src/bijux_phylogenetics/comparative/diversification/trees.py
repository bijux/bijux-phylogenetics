from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.diagnostics.validation import (
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    assess_tree_ultrametricity,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    DiversificationAnalysisError,
    NonUltrametricTreeError,
    UnrootedTreeError,
)

from .models import TimeTreeValidationReport


def node_depths(tree: PhyloTree) -> dict[str, float]:
    depths: dict[str, float] = {node_signature(tree.root): 0.0}

    def visit(node: TreeNode, depth: float) -> None:
        for child in node.children:
            next_depth = depth + float(child.branch_length or 0.0)
            depths[node_signature(child)] = next_depth
            if not child.is_leaf():
                visit(child, next_depth)

    visit(tree.root, 0.0)
    return depths


def root_age(tree: PhyloTree) -> float:
    distances = [
        distance for _tip, distance in tree.root_to_tip_pairs() if distance is not None
    ]
    if not distances:
        raise DiversificationAnalysisError(
            "diversification analysis requires complete root-to-tip distances"
        )
    return float(format(max(distances), ".15g"))


def descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(descendant_taxa(child))
    return sorted(taxa)


def is_strictly_bifurcating(tree: PhyloTree) -> bool:
    return all(
        len(node.children) == 2 for node in tree.iter_nodes() if not node.is_leaf()
    )


def node_age(tree: PhyloTree, depths: dict[str, float], node: TreeNode) -> float:
    return float(format(root_age(tree) - depths[node_signature(node)], ".15g"))


def find_smallest_covering_node(
    tree: PhyloTree, taxa: set[str]
) -> tuple[TreeNode, list[str]]:
    best_node = tree.root
    best_taxa = descendant_taxa(tree.root)
    for node in tree.iter_nodes():
        candidate_taxa = descendant_taxa(node)
        candidate_set = set(candidate_taxa)
        if taxa <= candidate_set and len(candidate_set) <= len(best_taxa):
            best_node = node
            best_taxa = candidate_taxa
    return best_node, best_taxa


def validate_time_tree_for_diversification(tree_path: Path) -> TimeTreeValidationReport:
    """Validate the rooted ultrametric time-tree contract required for diversification analysis."""
    validation = validate_tree_path(tree_path, require_rooted=True)
    if validation.branch_length_status != "complete":
        raise DiversificationAnalysisError(
            "diversification analysis requires complete branch lengths"
        )
    ultrametric_report = assess_tree_ultrametricity(tree_path)
    if ultrametric_report.rooted is not True:
        raise UnrootedTreeError(f"tree is not rooted: {tree_path}")
    if not ultrametric_report.ultrametric:
        raise NonUltrametricTreeError(
            f"tree is not ultrametric within ape-style diversification tolerance: {tree_path}",
            details={
                "tolerance": ultrametric_report.tolerance,
                "criterion_name": ultrametric_report.criterion_name,
                "criterion_value": ultrametric_report.criterion_value,
                "max_tip_depth_deviation": ultrametric_report.max_tip_depth_deviation,
                "offending_taxa": list(ultrametric_report.offending_taxa),
            },
        )
    return TimeTreeValidationReport(
        tree_path=tree_path,
        rooted=validation.rooted,
        ultrametric=True,
        branch_length_status=validation.branch_length_status,
        tip_count=validation.tip_count,
        root_age=float(format(ultrametric_report.root_age, ".15g")),
        warnings=list(validation.warnings),
    )


def inspect_diversification_time_tree(tree_path: Path) -> TimeTreeValidationReport:
    """Inspect time-tree readiness with explicit diversification semantics."""
    inspection = inspect_tree_path(tree_path)
    if inspection.branch_length_status != "complete":
        raise DiversificationAnalysisError(
            "diversification analysis requires complete branch lengths"
        )
    if not inspection.rooted:
        raise DiversificationAnalysisError(
            "diversification analysis requires a rooted tree"
        )
    ultrametric_report = assess_tree_ultrametricity(tree_path)
    if not ultrametric_report.ultrametric:
        raise DiversificationAnalysisError(
            "diversification analysis requires an ultrametric time tree"
        )
    return TimeTreeValidationReport(
        tree_path=tree_path,
        rooted=inspection.rooted,
        ultrametric=True,
        branch_length_status=inspection.branch_length_status,
        tip_count=inspection.tip_count,
        root_age=float(format(ultrametric_report.root_age, ".15g")),
        warnings=list(inspection.warnings),
    )
