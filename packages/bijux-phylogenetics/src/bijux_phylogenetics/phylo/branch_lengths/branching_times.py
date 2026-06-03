from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.phylo.branch_lengths.node_depths import (
    compute_tree_node_depths,
)
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    assess_tree_ultrametricity,
)
from bijux_phylogenetics.runtime.errors import (
    NonUltrametricTreeError,
    UnrootedTreeError,
)


@dataclass(slots=True)
class TreeBranchingTimeRow:
    """One ape-style internal-node branching-time row."""

    node_id: int
    node_kind: str
    node_label: str | None
    descendant_taxa: list[str]
    node_depth: float
    branching_time: float


@dataclass(slots=True)
class TreeBranchingTimeReport:
    """Deterministic ape-style branching times for one rooted ultrametric tree."""

    tree_path: Path
    tip_labels: list[str]
    rooted: bool | None
    tolerance: float
    tree_is_ultrametric: bool
    root_age: float
    internal_node_count: int
    zero_branch_length_count: int
    minimum_tip_depth: float
    maximum_tip_depth: float
    max_tip_depth_deviation: float
    rows: list[TreeBranchingTimeRow]


def compute_tree_branching_times(
    path: Path,
    *,
    tolerance: float = 1e-12,
) -> TreeBranchingTimeReport:
    """Compute ape-style branching times for one rooted ultrametric tree."""
    node_depth_report = compute_tree_node_depths(path)
    ultrametric_report = assess_tree_ultrametricity(path, tolerance=tolerance)
    if ultrametric_report.rooted is not True:
        raise UnrootedTreeError(
            "branching-time calculations require a rooted tree",
            code="tree_branching_times_require_rooted_tree",
            details={"rooted": ultrametric_report.rooted},
        )
    if not ultrametric_report.ultrametric:
        raise NonUltrametricTreeError(
            "branching-time calculations require an ultrametric tree",
            code="tree_branching_times_require_ultrametric_tree",
            details={
                "minimum_tip_depth": ultrametric_report.minimum_tip_depth,
                "maximum_tip_depth": ultrametric_report.maximum_tip_depth,
                "max_tip_depth_deviation": ultrametric_report.max_tip_depth_deviation,
                "offending_taxa": list(ultrametric_report.offending_taxa),
                "tolerance": ultrametric_report.tolerance,
            },
        )
    root_age = ultrametric_report.root_age
    rows = [
        TreeBranchingTimeRow(
            node_id=row.node_id,
            node_kind=row.node_kind,
            node_label=row.node_label,
            descendant_taxa=row.descendant_taxa,
            node_depth=row.branch_length_depth,
            branching_time=root_age - row.branch_length_depth,
        )
        for row in node_depth_report.rows
        if row.node_kind != "tip"
    ]
    return TreeBranchingTimeReport(
        tree_path=path,
        tip_labels=node_depth_report.tip_labels,
        rooted=ultrametric_report.rooted,
        tolerance=tolerance,
        tree_is_ultrametric=True,
        root_age=root_age,
        internal_node_count=node_depth_report.internal_node_count,
        zero_branch_length_count=node_depth_report.zero_branch_length_count,
        minimum_tip_depth=ultrametric_report.minimum_tip_depth,
        maximum_tip_depth=ultrametric_report.maximum_tip_depth,
        max_tip_depth_deviation=ultrametric_report.max_tip_depth_deviation,
        rows=rows,
    )


def write_tree_branching_time_table(
    path: Path,
    report: TreeBranchingTimeReport,
) -> Path:
    """Write one deterministic ape-style branching-time ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        (
            "node_id\tnode_kind\tnode_label\tdescendant_taxa\tnode_depth\t"
            "branching_time\trooted\troot_age\ttree_is_ultrametric\t"
            "internal_node_count\tzero_branch_length_count\tminimum_tip_depth\t"
            "maximum_tip_depth\tmax_tip_depth_deviation\ttolerance"
        )
    ]
    for row in report.rows:
        lines.append(
            "\t".join(
                [
                    str(row.node_id),
                    row.node_kind,
                    row.node_label or "",
                    "|".join(row.descendant_taxa),
                    format(row.node_depth, ".15g"),
                    format(row.branching_time, ".15g"),
                    str(report.rooted),
                    format(report.root_age, ".15g"),
                    str(report.tree_is_ultrametric),
                    str(report.internal_node_count),
                    str(report.zero_branch_length_count),
                    format(report.minimum_tip_depth, ".15g"),
                    format(report.maximum_tip_depth, ".15g"),
                    format(report.max_tip_depth_deviation, ".15g"),
                    format(report.tolerance, ".15g"),
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
