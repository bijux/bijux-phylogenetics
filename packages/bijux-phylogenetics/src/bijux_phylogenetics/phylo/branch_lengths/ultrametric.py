from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import sys

from bijux_phylogenetics.phylo.branch_lengths.node_depths import (
    TreeNodeDepthReport,
    compute_tree_node_depths,
)

APE_ULTRAMETRIC_TOLERANCE = sys.float_info.epsilon**0.5


@dataclass(slots=True)
class TreeUltrametricTipRow:
    """One tip-level root-to-tip depth diagnostic row."""

    node_id: int
    tip_label: str
    root_to_tip_depth: float
    deviation_from_mean_depth: float
    deviation_from_min_depth: float
    deviation_from_max_depth: float
    is_offending_taxon: bool


@dataclass(slots=True)
class TreeUltrametricReport:
    """One governed ape-style ultrametricity assessment."""

    tree_path: Path
    tip_labels: list[str]
    rooted: bool | None
    tolerance: float
    option: int
    criterion_name: str
    criterion_value: float
    ultrametric: bool
    minimum_tip_depth: float
    maximum_tip_depth: float
    mean_tip_depth: float
    max_tip_depth_deviation: float
    root_age: float
    offending_taxa: list[str]
    rows: list[TreeUltrametricTipRow]


@dataclass(slots=True)
class TipDepthUltrametricSummary:
    """Shared ultrametricity summary derived from labeled tip depths."""

    tolerance: float
    option: int
    criterion_name: str
    criterion_value: float
    ultrametric: bool
    minimum_tip_depth: float
    maximum_tip_depth: float
    mean_tip_depth: float
    max_tip_depth_deviation: float
    root_age: float
    offending_taxa: list[str]


def assess_tree_ultrametricity(
    path: Path,
    *,
    tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
    option: int = 1,
) -> TreeUltrametricReport:
    """Assess ape-style ultrametricity from root-to-tip branch-length depths."""
    node_depth_report = compute_tree_node_depths(path)
    return _summarize_tree_ultrametricity(
        node_depth_report,
        tree_path=path,
        tolerance=tolerance,
        option=option,
    )


def summarize_ultrametric_tip_depths(
    tip_depth_by_label: dict[str, float],
    *,
    tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
    option: int = 1,
) -> TipDepthUltrametricSummary:
    """Summarize ultrametricity from one explicit tip-depth mapping."""
    if option not in {1, 2}:
        raise ValueError("ultrametric option must be 1 or 2")
    tip_depths = list(tip_depth_by_label.values())
    minimum_tip_depth = min(tip_depths)
    maximum_tip_depth = max(tip_depths)
    mean_tip_depth = sum(tip_depths) / len(tip_depths)
    max_tip_depth_deviation = maximum_tip_depth - minimum_tip_depth
    criterion_name = "scaled-range" if option == 1 else "variance"
    criterion_value = _criterion_value(
        tip_depths,
        maximum_tip_depth=maximum_tip_depth,
        max_tip_depth_deviation=max_tip_depth_deviation,
        option=option,
    )
    offending_taxa = _offending_taxa_by_depth(
        tip_depth_by_label,
        minimum_tip_depth=minimum_tip_depth,
        maximum_tip_depth=maximum_tip_depth,
    )
    return TipDepthUltrametricSummary(
        tolerance=tolerance,
        option=option,
        criterion_name=criterion_name,
        criterion_value=criterion_value,
        ultrametric=criterion_value <= tolerance,
        minimum_tip_depth=minimum_tip_depth,
        maximum_tip_depth=maximum_tip_depth,
        mean_tip_depth=mean_tip_depth,
        max_tip_depth_deviation=max_tip_depth_deviation,
        root_age=maximum_tip_depth,
        offending_taxa=offending_taxa,
    )


def write_tree_ultrametric_table(
    path: Path,
    report: TreeUltrametricReport,
) -> Path:
    """Write one deterministic ultrametricity diagnostic ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        (
            "node_id\ttip_label\troot_to_tip_depth\tdeviation_from_mean_depth\t"
            "deviation_from_min_depth\tdeviation_from_max_depth\t"
            "is_offending_taxon\trooted\tultrametric\tcriterion_name\t"
            "criterion_value\ttolerance\tminimum_tip_depth\tmaximum_tip_depth\t"
            "mean_tip_depth\tmax_tip_depth_deviation\troot_age"
        )
    ]
    for row in report.rows:
        lines.append(
            "\t".join(
                [
                    str(row.node_id),
                    row.tip_label,
                    format(row.root_to_tip_depth, ".15g"),
                    format(row.deviation_from_mean_depth, ".15g"),
                    format(row.deviation_from_min_depth, ".15g"),
                    format(row.deviation_from_max_depth, ".15g"),
                    str(row.is_offending_taxon).lower(),
                    str(report.rooted),
                    str(report.ultrametric),
                    report.criterion_name,
                    format(report.criterion_value, ".15g"),
                    format(report.tolerance, ".15g"),
                    format(report.minimum_tip_depth, ".15g"),
                    format(report.maximum_tip_depth, ".15g"),
                    format(report.mean_tip_depth, ".15g"),
                    format(report.max_tip_depth_deviation, ".15g"),
                    format(report.root_age, ".15g"),
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _summarize_tree_ultrametricity(
    node_depth_report: TreeNodeDepthReport,
    *,
    tree_path: Path,
    tolerance: float,
    option: int,
) -> TreeUltrametricReport:
    tip_rows = [row for row in node_depth_report.rows if row.node_kind == "tip"]
    summary = summarize_ultrametric_tip_depths(
        {row.node_label or "": row.branch_length_depth for row in tip_rows},
        tolerance=tolerance,
        option=option,
    )
    return TreeUltrametricReport(
        tree_path=tree_path,
        tip_labels=node_depth_report.tip_labels,
        rooted=node_depth_report.rooted,
        tolerance=summary.tolerance,
        option=summary.option,
        criterion_name=summary.criterion_name,
        criterion_value=summary.criterion_value,
        ultrametric=summary.ultrametric,
        minimum_tip_depth=summary.minimum_tip_depth,
        maximum_tip_depth=summary.maximum_tip_depth,
        mean_tip_depth=summary.mean_tip_depth,
        max_tip_depth_deviation=summary.max_tip_depth_deviation,
        root_age=summary.root_age,
        offending_taxa=summary.offending_taxa,
        rows=[
            TreeUltrametricTipRow(
                node_id=row.node_id,
                tip_label=row.node_label or "",
                root_to_tip_depth=row.branch_length_depth,
                deviation_from_mean_depth=abs(
                    row.branch_length_depth - summary.mean_tip_depth
                ),
                deviation_from_min_depth=(
                    row.branch_length_depth - summary.minimum_tip_depth
                ),
                deviation_from_max_depth=(
                    summary.maximum_tip_depth - row.branch_length_depth
                ),
                is_offending_taxon=(row.node_label or "")
                in set(summary.offending_taxa),
            )
            for row in tip_rows
        ],
    )


def _criterion_value(
    tip_depths: list[float],
    *,
    maximum_tip_depth: float,
    max_tip_depth_deviation: float,
    option: int,
) -> float:
    if option == 1:
        if math.isclose(maximum_tip_depth, 0.0, abs_tol=1e-15):
            if math.isclose(max_tip_depth_deviation, 0.0, abs_tol=1e-15):
                return 0.0
            return math.inf
        return max_tip_depth_deviation / maximum_tip_depth
    if len(tip_depths) <= 1:
        return 0.0
    mean_tip_depth = sum(tip_depths) / len(tip_depths)
    squared_error = sum((depth - mean_tip_depth) ** 2 for depth in tip_depths)
    return squared_error / (len(tip_depths) - 1)


def _offending_taxa_by_depth(
    tip_depth_by_label: dict[str, float],
    *,
    minimum_tip_depth: float,
    maximum_tip_depth: float,
) -> list[str]:
    if math.isclose(minimum_tip_depth, maximum_tip_depth, abs_tol=1e-12):
        return []
    offending: list[str] = []
    for label, depth in tip_depth_by_label.items():
        if label == "":
            continue
        if math.isclose(depth, minimum_tip_depth, abs_tol=1e-12):
            offending.append(label)
            continue
        if math.isclose(depth, maximum_tip_depth, abs_tol=1e-12):
            offending.append(label)
    return sorted(set(offending))
