from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeReadinessReport,
    node_signature,
    summarize_numeric_trait_readiness,
    tip_root_depths,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


@dataclass(slots=True)
class TraitRateThroughTimeExclusion:
    """One taxon excluded before rate-through-time trait analysis."""

    taxon: str
    reason: str


@dataclass(slots=True)
class TraitRateThroughTimeIntervalRow:
    """One time interval with aggregated branchwise trait-rate evidence."""

    interval_index: int
    start_depth: float
    end_depth: float
    midpoint_depth: float
    overlapping_branch_count: int
    segment_count: int
    total_segment_length: float
    total_squared_change: float
    estimated_rate: float | None
    share_of_total_squared_change: float | None


@dataclass(slots=True)
class TraitRateThroughTimeSummaryReport:
    """Reviewer-facing rate-through-time summary for one continuous trait."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    analyzed_taxon_count: int
    excluded_taxa: list[TraitRateThroughTimeExclusion]
    interval_count: int
    nonempty_interval_count: int
    tree_depth: float
    interval_rows: list[TraitRateThroughTimeIntervalRow]
    ancestral_model: str
    earliest_interval_rate: float | None
    latest_interval_rate: float | None
    latest_to_earliest_rate_ratio: float | None
    weighted_rate_slope: float | None
    normalized_rate_slope: float | None
    trend_direction: str
    peak_interval_index: int | None
    trough_interval_index: int | None
    assumptions: list[str]
    warnings: list[str]
    readiness: ComparativeReadinessReport


def summarize_trait_rate_through_time(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    interval_count: int = 5,
) -> TraitRateThroughTimeSummaryReport:
    """Summarize how reconstructed branchwise trait change varies across tree depth."""
    if interval_count < 2:
        raise ComparativeMethodError(
            "rate-through-time trait analysis requires at least two intervals"
        )
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    if not readiness.rooted:
        raise ComparativeMethodError(
            "rate-through-time trait analysis requires a rooted tree"
        )
    if not readiness.binary:
        raise ComparativeMethodError(
            "rate-through-time trait analysis requires a strictly binary tree"
        )
    if not readiness.complete_branch_lengths:
        raise ComparativeMethodError(
            "rate-through-time trait analysis requires complete branch lengths"
        )
    if len(readiness.analysis_taxa) < 3:
        raise ComparativeMethodError(
            "rate-through-time trait analysis requires at least three taxa with usable numeric trait values"
        )
    pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, readiness.analysis_taxa)
    if any(
        (node.branch_length or 0.0) <= 0.0
        for node in pruned_tree.iter_nodes()
        if node is not pruned_tree.root
    ):
        raise ComparativeMethodError(
            "rate-through-time trait analysis requires strictly positive branch lengths"
        )
    tree_depth = max(tip_root_depths(pruned_tree, pruned_tree.tip_names).values())
    if tree_depth <= 0.0:
        raise ComparativeMethodError(
            "rate-through-time trait analysis requires a positive total tree depth"
        )
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-rate-through-time-"
    ) as tmp_dir:
        filtered_table_path = Path(tmp_dir) / "rate-through-time-traits.tsv"
        write_taxon_rows(
            filtered_table_path,
            columns=[readiness.taxon_column, trait],
            rows=_build_filtered_trait_rows(
                traits_path=traits_path,
                taxon_column=readiness.taxon_column,
                trait=trait,
                analyzed_taxa=readiness.analysis_taxa,
            ),
        )
        reconstruction = reconstruct_continuous_ancestral_states(
            tree_path,
            filtered_table_path,
            trait=trait,
            taxon_column=readiness.taxon_column,
            model="brownian",
        )
    estimate_lookup = {
        estimate.node: estimate.estimate for estimate in reconstruction.estimates
    }
    interval_rows = _build_interval_rows(
        pruned_tree=pruned_tree,
        estimate_lookup=estimate_lookup,
        tree_depth=tree_depth,
        interval_count=interval_count,
    )
    nonempty_rows = [row for row in interval_rows if row.estimated_rate is not None]
    earliest_rate = nonempty_rows[0].estimated_rate if nonempty_rows else None
    latest_rate = nonempty_rows[-1].estimated_rate if nonempty_rows else None
    rate_ratio = _rate_ratio(latest_rate, earliest_rate)
    slope = _weighted_rate_slope(nonempty_rows)
    normalized_slope = _normalized_rate_slope(
        slope=slope,
        tree_depth=tree_depth,
        rows=nonempty_rows,
    )
    trend_direction = _classify_trend_direction(
        earliest_rate=earliest_rate,
        latest_rate=latest_rate,
        rate_ratio=rate_ratio,
        normalized_slope=normalized_slope,
        nonempty_interval_count=len(nonempty_rows),
    )
    peak_interval_index = (
        max(
            nonempty_rows, key=lambda row: row.estimated_rate or -math.inf
        ).interval_index
        if nonempty_rows
        else None
    )
    trough_interval_index = (
        min(
            nonempty_rows, key=lambda row: row.estimated_rate or math.inf
        ).interval_index
        if nonempty_rows
        else None
    )
    warnings = list(dict.fromkeys([*readiness.warnings, *reconstruction.warnings]))
    if any(row.estimated_rate is None for row in interval_rows):
        warnings.append(
            "one or more rate-through-time intervals contain no branch overlap on the analyzed tree"
        )
    assumptions = [
        "Branchwise trait-rate evidence is derived from squared change between Brownian ancestral estimates and observed or reconstructed descendant values.",
        "When a branch spans more than one interval, its squared change is allocated across intervals in proportion to overlapping branch length.",
        "Trend classification is reviewer-facing and based on weighted interval-rate slope plus the latest-to-earliest rate ratio.",
    ]
    return TraitRateThroughTimeSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=reconstruction.taxon_column,
        trait=trait,
        tree_taxon_count=readiness.tree_taxa,
        analyzed_taxa=list(readiness.analysis_taxa),
        analyzed_taxon_count=len(readiness.analysis_taxa),
        excluded_taxa=_build_excluded_taxa(readiness),
        interval_count=interval_count,
        nonempty_interval_count=len(nonempty_rows),
        tree_depth=tree_depth,
        interval_rows=interval_rows,
        ancestral_model="brownian",
        earliest_interval_rate=earliest_rate,
        latest_interval_rate=latest_rate,
        latest_to_earliest_rate_ratio=rate_ratio,
        weighted_rate_slope=slope,
        normalized_rate_slope=normalized_slope,
        trend_direction=trend_direction,
        peak_interval_index=peak_interval_index,
        trough_interval_index=trough_interval_index,
        assumptions=assumptions,
        warnings=warnings,
        readiness=readiness,
    )


def write_trait_rate_through_time_summary_table(
    path: Path,
    report: TraitRateThroughTimeSummaryReport,
) -> Path:
    """Write one summary ledger for rate-through-time trait analysis."""
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "interval_count",
            "nonempty_interval_count",
            "tree_depth",
            "ancestral_model",
            "earliest_interval_rate",
            "latest_interval_rate",
            "latest_to_earliest_rate_ratio",
            "weighted_rate_slope",
            "normalized_rate_slope",
            "trend_direction",
            "peak_interval_index",
            "trough_interval_index",
            "warning_count",
        ],
        rows=[
            {
                "trait": report.trait,
                "taxon_column": report.taxon_column,
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": report.analyzed_taxon_count,
                "excluded_taxon_count": len(report.excluded_taxa),
                "interval_count": report.interval_count,
                "nonempty_interval_count": report.nonempty_interval_count,
                "tree_depth": format(report.tree_depth, ".15g"),
                "ancestral_model": report.ancestral_model,
                "earliest_interval_rate": _format_optional(
                    report.earliest_interval_rate
                ),
                "latest_interval_rate": _format_optional(report.latest_interval_rate),
                "latest_to_earliest_rate_ratio": _format_optional(
                    report.latest_to_earliest_rate_ratio
                ),
                "weighted_rate_slope": _format_optional(report.weighted_rate_slope),
                "normalized_rate_slope": _format_optional(report.normalized_rate_slope),
                "trend_direction": report.trend_direction,
                "peak_interval_index": (
                    ""
                    if report.peak_interval_index is None
                    else report.peak_interval_index
                ),
                "trough_interval_index": (
                    ""
                    if report.trough_interval_index is None
                    else report.trough_interval_index
                ),
                "warning_count": len(report.warnings),
            }
        ],
    )


def write_trait_rate_through_time_interval_table(
    path: Path,
    report: TraitRateThroughTimeSummaryReport,
) -> Path:
    """Write one interval ledger for rate-through-time trait analysis."""
    return write_taxon_rows(
        path,
        columns=[
            "interval_index",
            "start_depth",
            "end_depth",
            "midpoint_depth",
            "overlapping_branch_count",
            "segment_count",
            "total_segment_length",
            "total_squared_change",
            "estimated_rate",
            "share_of_total_squared_change",
        ],
        rows=[
            {
                "interval_index": row.interval_index,
                "start_depth": format(row.start_depth, ".15g"),
                "end_depth": format(row.end_depth, ".15g"),
                "midpoint_depth": format(row.midpoint_depth, ".15g"),
                "overlapping_branch_count": row.overlapping_branch_count,
                "segment_count": row.segment_count,
                "total_segment_length": format(row.total_segment_length, ".15g"),
                "total_squared_change": format(row.total_squared_change, ".15g"),
                "estimated_rate": _format_optional(row.estimated_rate),
                "share_of_total_squared_change": _format_optional(
                    row.share_of_total_squared_change
                ),
            }
            for row in report.interval_rows
        ],
    )


def write_trait_rate_through_time_exclusion_table(
    path: Path,
    report: TraitRateThroughTimeSummaryReport,
) -> Path:
    """Write one excluded-taxon ledger for rate-through-time trait analysis."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {"taxon": row.taxon, "reason": row.reason} for row in report.excluded_taxa
        ],
    )


@dataclass(slots=True)
class _IntervalAccumulator:
    interval_index: int
    start_depth: float
    end_depth: float
    midpoint_depth: float
    overlapping_branch_ids: set[str]
    segment_count: int
    total_segment_length: float
    total_squared_change: float


def _build_interval_rows(
    *,
    pruned_tree: PhyloTree,
    estimate_lookup: dict[str, float],
    tree_depth: float,
    interval_count: int,
) -> list[TraitRateThroughTimeIntervalRow]:
    width = tree_depth / interval_count
    intervals = [
        _IntervalAccumulator(
            interval_index=index + 1,
            start_depth=index * width,
            end_depth=tree_depth
            if index == interval_count - 1
            else (index + 1) * width,
            midpoint_depth=(
                (index * width)
                + (tree_depth if index == interval_count - 1 else (index + 1) * width)
            )
            / 2.0,
            overlapping_branch_ids=set(),
            segment_count=0,
            total_segment_length=0.0,
            total_squared_change=0.0,
        )
        for index in range(interval_count)
    ]
    root_signature = node_signature(pruned_tree.root)
    root_estimate = estimate_lookup[root_signature]

    def visit(node: TreeNode, *, parent_depth: float, parent_estimate: float) -> None:
        branch_length = float(node.branch_length or 0.0)
        child_depth = parent_depth + branch_length
        branch_id = node_signature(node)
        child_estimate = estimate_lookup[branch_id]
        squared_change = (child_estimate - parent_estimate) ** 2
        for interval in intervals:
            overlap = min(child_depth, interval.end_depth) - max(
                parent_depth, interval.start_depth
            )
            if overlap <= 0.0:
                continue
            interval.overlapping_branch_ids.add(branch_id)
            interval.segment_count += 1
            interval.total_segment_length += overlap
            interval.total_squared_change += squared_change * (overlap / branch_length)
        for child in node.children:
            visit(
                child,
                parent_depth=child_depth,
                parent_estimate=child_estimate,
            )

    for child in pruned_tree.root.children:
        visit(child, parent_depth=0.0, parent_estimate=root_estimate)

    total_squared_change = sum(interval.total_squared_change for interval in intervals)
    return [
        TraitRateThroughTimeIntervalRow(
            interval_index=interval.interval_index,
            start_depth=interval.start_depth,
            end_depth=interval.end_depth,
            midpoint_depth=interval.midpoint_depth,
            overlapping_branch_count=len(interval.overlapping_branch_ids),
            segment_count=interval.segment_count,
            total_segment_length=interval.total_segment_length,
            total_squared_change=interval.total_squared_change,
            estimated_rate=(
                None
                if interval.total_segment_length <= 0.0
                else interval.total_squared_change / interval.total_segment_length
            ),
            share_of_total_squared_change=(
                None
                if total_squared_change <= 0.0
                else interval.total_squared_change / total_squared_change
            ),
        )
        for interval in intervals
    ]


def _build_excluded_taxa(
    readiness: ComparativeReadinessReport,
) -> list[TraitRateThroughTimeExclusion]:
    rows: list[TraitRateThroughTimeExclusion] = []
    for taxon in readiness.missing_from_traits:
        rows.append(
            TraitRateThroughTimeExclusion(
                taxon=taxon,
                reason="missing_from_trait_table",
            )
        )
    for taxon in readiness.pruned_missing_value_taxa:
        rows.append(
            TraitRateThroughTimeExclusion(
                taxon=taxon,
                reason="missing_trait_value",
            )
        )
    for taxon in readiness.pruned_non_numeric_taxa:
        rows.append(
            TraitRateThroughTimeExclusion(
                taxon=taxon,
                reason="non_numeric_trait_value",
            )
        )
    for taxon in readiness.extra_trait_taxa:
        rows.append(
            TraitRateThroughTimeExclusion(
                taxon=taxon,
                reason="absent_from_tree",
            )
        )
    return rows


def _build_filtered_trait_rows(
    *,
    traits_path: Path,
    taxon_column: str,
    trait: str,
    analyzed_taxa: list[str],
) -> list[dict[str, str]]:
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    wanted = set(analyzed_taxa)
    return [
        {
            taxon_column: row[taxon_column],
            trait: row[trait],
        }
        for row in table.rows
        if row[taxon_column] in wanted
    ]


def _rate_ratio(
    latest_rate: float | None,
    earliest_rate: float | None,
) -> float | None:
    if latest_rate is None or earliest_rate is None or earliest_rate <= 0.0:
        return None
    return latest_rate / earliest_rate


def _weighted_rate_slope(
    rows: list[TraitRateThroughTimeIntervalRow],
) -> float | None:
    if len(rows) < 2:
        return None
    weights = [row.total_segment_length for row in rows]
    total_weight = sum(weights)
    if total_weight <= 0.0:
        return None
    mean_x = (
        sum(
            weight * row.midpoint_depth
            for weight, row in zip(weights, rows, strict=True)
        )
        / total_weight
    )
    mean_y = (
        sum(
            weight * (row.estimated_rate or 0.0)
            for weight, row in zip(weights, rows, strict=True)
        )
        / total_weight
    )
    denominator = sum(
        weight * (row.midpoint_depth - mean_x) ** 2
        for weight, row in zip(weights, rows, strict=True)
    )
    if denominator <= 1e-12:
        return 0.0
    numerator = sum(
        weight * (row.midpoint_depth - mean_x) * ((row.estimated_rate or 0.0) - mean_y)
        for weight, row in zip(weights, rows, strict=True)
    )
    return numerator / denominator


def _normalized_rate_slope(
    *,
    slope: float | None,
    tree_depth: float,
    rows: list[TraitRateThroughTimeIntervalRow],
) -> float | None:
    if slope is None or not rows:
        return None
    total_length = sum(row.total_segment_length for row in rows)
    if total_length <= 0.0:
        return None
    mean_rate = (
        sum((row.estimated_rate or 0.0) * row.total_segment_length for row in rows)
        / total_length
    )
    if mean_rate <= 0.0:
        return None
    return (slope * tree_depth) / mean_rate


def _classify_trend_direction(
    *,
    earliest_rate: float | None,
    latest_rate: float | None,
    rate_ratio: float | None,
    normalized_slope: float | None,
    nonempty_interval_count: int,
) -> str:
    if nonempty_interval_count < 2 or normalized_slope is None or rate_ratio is None:
        return "insufficient_data"
    if (
        earliest_rate is None
        or latest_rate is None
        or earliest_rate <= 0.0
        or latest_rate <= 0.0
    ):
        return "insufficient_data"
    if normalized_slope <= -0.2 and rate_ratio <= 0.85:
        return "slowdown"
    if normalized_slope >= 0.2 and rate_ratio >= 1.15:
        return "acceleration"
    return "stable"


def _format_optional(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")
