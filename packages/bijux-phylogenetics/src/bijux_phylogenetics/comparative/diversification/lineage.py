from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree

from .models import LineageThroughTimePoint, LineageThroughTimeReport
from .trees import node_depths, validate_time_tree_for_diversification


def compute_lineage_through_time_curve(tree_path: Path) -> LineageThroughTimeReport:
    """Compute a deterministic lineage-through-time summary for one valid time tree."""
    validation = validate_time_tree_for_diversification(tree_path)
    tree = load_tree(tree_path)
    depths = node_depths(tree)
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


def write_lineage_through_time_table(
    path: Path, report: LineageThroughTimeReport
) -> Path:
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
