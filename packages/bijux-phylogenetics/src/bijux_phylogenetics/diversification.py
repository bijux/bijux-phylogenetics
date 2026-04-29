from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.diagnostics.validation import inspect_tree_path, validate_tree_path
from bijux_phylogenetics.errors import DiversificationAnalysisError
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class TimeTreeValidationReport:
    tree_path: Path
    rooted: bool
    ultrametric: bool
    branch_length_status: str
    tip_count: int
    root_age: float
    warnings: list[str]


@dataclass(slots=True)
class LineageThroughTimePoint:
    node: str
    time_before_present: float
    lineage_count: int
    event: str


@dataclass(slots=True)
class LineageThroughTimeReport:
    tree_path: Path
    tip_count: int
    root_age: float
    points: list[LineageThroughTimePoint]


def _node_depths(tree: PhyloTree) -> dict[str, float]:
    depths: dict[str, float] = {node_signature(tree.root): 0.0}

    def visit(node: TreeNode, depth: float) -> None:
        for child in node.children:
            next_depth = depth + float(child.branch_length or 0.0)
            depths[node_signature(child)] = next_depth
            if not child.is_leaf():
                visit(child, next_depth)

    visit(tree.root, 0.0)
    return depths


def _root_age(tree: PhyloTree) -> float:
    distances = [distance for _tip, distance in tree.root_to_tip_pairs() if distance is not None]
    if not distances:
        raise DiversificationAnalysisError("diversification analysis requires complete root-to-tip distances")
    return float(format(max(distances), ".15g"))


def validate_time_tree_for_diversification(tree_path: Path) -> TimeTreeValidationReport:
    """Validate the rooted ultrametric time-tree contract required for diversification analysis."""
    validation = validate_tree_path(
        tree_path,
        require_rooted=True,
        require_ultrametric=True,
    )
    if validation.branch_length_status != "complete":
        raise DiversificationAnalysisError("diversification analysis requires complete branch lengths")
    tree = load_tree(tree_path)
    root_age = _root_age(tree)
    return TimeTreeValidationReport(
        tree_path=tree_path,
        rooted=validation.rooted,
        ultrametric=validation.ultrametric is True,
        branch_length_status=validation.branch_length_status,
        tip_count=validation.tip_count,
        root_age=root_age,
        warnings=list(validation.warnings),
    )


def compute_lineage_through_time_curve(tree_path: Path) -> LineageThroughTimeReport:
    """Compute a deterministic lineage-through-time summary for one valid time tree."""
    validation = validate_time_tree_for_diversification(tree_path)
    tree = load_tree(tree_path)
    depths = _node_depths(tree)
    root_age = validation.root_age
    events = [
        (
            float(format(root_age - depths[node_signature(node)], ".15g")),
            node_signature(node),
            max(len(node.children) - 1, 0),
        )
        for node in tree.iter_nodes()
        if node is not tree.root and not node.is_leaf()
    ]
    events.sort(key=lambda item: (-item[0], item[1]))

    lineage_count = max(len(tree.root.children), 1)
    points = [
        LineageThroughTimePoint(
            node=node_signature(tree.root),
            time_before_present=root_age,
            lineage_count=lineage_count,
            event="root",
        )
    ]
    for event_age, signature, gained_lineages in events:
        lineage_count += gained_lineages
        points.append(
            LineageThroughTimePoint(
                node=signature,
                time_before_present=event_age,
                lineage_count=lineage_count,
                event="speciation",
            )
        )
    if points[-1].time_before_present != 0.0:
        points.append(
            LineageThroughTimePoint(
                node="present",
                time_before_present=0.0,
                lineage_count=tree.tip_count,
                event="present",
            )
        )
    return LineageThroughTimeReport(
        tree_path=tree_path,
        tip_count=tree.tip_count,
        root_age=root_age,
        points=points,
    )


def write_lineage_through_time_table(path: Path, report: LineageThroughTimeReport) -> Path:
    """Export a lineage-through-time curve as a deterministic table."""
    rows = [
        {
            "node": point.node,
            "time_before_present": format(point.time_before_present, ".15g"),
            "lineage_count": str(point.lineage_count),
            "event": point.event,
        }
        for point in report.points
    ]
    return write_taxon_rows(
        path,
        columns=["node", "time_before_present", "lineage_count", "event"],
        rows=rows,
    )


def inspect_diversification_time_tree(tree_path: Path) -> TimeTreeValidationReport:
    """Inspect time-tree readiness with explicit diversification semantics."""
    inspection = inspect_tree_path(tree_path)
    if inspection.branch_length_status != "complete":
        raise DiversificationAnalysisError("diversification analysis requires complete branch lengths")
    if not inspection.rooted:
        raise DiversificationAnalysisError("diversification analysis requires a rooted tree")
    if inspection.is_ultrametric is not True:
        raise DiversificationAnalysisError("diversification analysis requires an ultrametric time tree")
    tree = load_tree(tree_path)
    return TimeTreeValidationReport(
        tree_path=tree_path,
        rooted=inspection.rooted,
        ultrametric=inspection.is_ultrametric is True,
        branch_length_status=inspection.branch_length_status,
        tip_count=inspection.tip_count,
        root_age=_root_age(tree),
        warnings=list(inspection.warnings),
    )
