from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_signature, stable_value
from bijux_phylogenetics.biogeography.state_models import (
    GeographicExcludedTaxonRow,
    GeographicStateModelReport,
    summarize_geographic_state_model,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError


@dataclass(frozen=True, slots=True)
class TimeBinDefinition:
    """One explicit root-depth interval for time-stratified geography."""

    label: str
    start_depth: float
    end_depth: float


@dataclass(frozen=True, slots=True)
class TimeStratifiedTransitionMatrixRow:
    """One source-target transition rate inside one time bin."""

    time_bin_label: str
    start_depth: float
    end_depth: float
    source_region: str
    target_region: str
    source_exposure_length: float
    allocated_transition_weight: float
    time_stratified_rate: float
    global_rate: float


@dataclass(frozen=True, slots=True)
class TimeStratifiedBranchRow:
    """One branch overlap with one time bin."""

    time_bin_label: str
    start_depth: float
    end_depth: float
    parent_node: str
    child_node: str
    parent_depth: float
    child_depth: float
    source_region: str
    target_region: str
    changed: bool
    overlap_length: float
    allocated_transition_weight: float
    support: float
    strongly_supported: bool


@dataclass(frozen=True, slots=True)
class TimeStratifiedTransitionSummary:
    """One reviewer-facing summary row for time-stratified geography."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    likelihood_method: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    tree_depth: float
    time_bin_count: int
    matrix_row_count: int
    changed_branch_count: int
    allocated_transition_weight_total: float
    warning_count: int


@dataclass(slots=True)
class TimeStratifiedTransitionReport:
    """Time-stratified geographic transition review surface."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    likelihood_method: str
    time_bins: list[TimeBinDefinition]
    summary: TimeStratifiedTransitionSummary
    matrix_rows: list[TimeStratifiedTransitionMatrixRow]
    branch_rows: list[TimeStratifiedBranchRow]
    exclusion_rows: list[GeographicExcludedTaxonRow]
    warnings: list[str]


def summarize_time_stratified_geographic_transitions(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    time_bins: list[TimeBinDefinition],
    taxon_column: str | None = None,
    model: str = "er",
    allowed_regions: list[str] | None = None,
) -> TimeStratifiedTransitionReport:
    """Estimate interval-specific geographic transitions from one rooted tree."""
    resolved_bins = _normalize_time_bins(time_bins)
    base_report = summarize_geographic_state_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_regions=allowed_regions,
    )
    tree = load_tree(tree_path)
    depth_by_node = _node_depths(tree)
    tree_depth = _tree_depth(tree)
    _validate_time_stratified_tree(tree, tree_depth)
    _validate_time_bins(resolved_bins, tree_depth)
    branch_rows = _build_branch_rows(base_report, resolved_bins, depth_by_node)
    matrix_rows = _build_matrix_rows(base_report, resolved_bins, branch_rows)
    warnings = list(base_report.warnings)
    coverage_warning = _build_time_bin_coverage_warning(resolved_bins, tree_depth)
    if coverage_warning is not None:
        warnings.append(coverage_warning)
    summary = TimeStratifiedTransitionSummary(
        trait=base_report.trait,
        taxon_column=base_report.taxon_column,
        model=base_report.model,
        internal_model=base_report.internal_model,
        likelihood_method=base_report.likelihood_method,
        analyzed_taxon_count=base_report.summary.analyzed_taxon_count,
        excluded_taxon_count=base_report.summary.excluded_taxon_count,
        tree_depth=tree_depth,
        time_bin_count=len(resolved_bins),
        matrix_row_count=len(matrix_rows),
        changed_branch_count=sum(
            row.changed for row in base_report.transition_event_rows
        ),
        allocated_transition_weight_total=stable_value(
            sum(row.allocated_transition_weight for row in branch_rows)
        ),
        warning_count=len(warnings),
    )
    return TimeStratifiedTransitionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=base_report.trait,
        taxon_column=base_report.taxon_column,
        model=base_report.model,
        internal_model=base_report.internal_model,
        likelihood_method=base_report.likelihood_method,
        time_bins=resolved_bins,
        summary=summary,
        matrix_rows=matrix_rows,
        branch_rows=branch_rows,
        exclusion_rows=list(base_report.exclusion_rows),
        warnings=warnings,
    )


def write_time_stratified_transition_summary_table(
    path: Path,
    report: TimeStratifiedTransitionReport,
) -> Path:
    """Write one summary ledger for time-stratified geographic transitions."""
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "internal_model",
            "likelihood_method",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "tree_depth",
            "time_bin_count",
            "matrix_row_count",
            "changed_branch_count",
            "allocated_transition_weight_total",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "internal_model": summary.internal_model,
                "likelihood_method": summary.likelihood_method,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "tree_depth": str(summary.tree_depth),
                "time_bin_count": str(summary.time_bin_count),
                "matrix_row_count": str(summary.matrix_row_count),
                "changed_branch_count": str(summary.changed_branch_count),
                "allocated_transition_weight_total": str(
                    summary.allocated_transition_weight_total
                ),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_time_stratified_transition_matrix_table(
    path: Path,
    report: TimeStratifiedTransitionReport,
) -> Path:
    """Write one time-specific transition matrix ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "time_bin_label",
            "start_depth",
            "end_depth",
            "source_region",
            "target_region",
            "source_exposure_length",
            "allocated_transition_weight",
            "time_stratified_rate",
            "global_rate",
        ],
        rows=[
            {
                "time_bin_label": row.time_bin_label,
                "start_depth": str(row.start_depth),
                "end_depth": str(row.end_depth),
                "source_region": row.source_region,
                "target_region": row.target_region,
                "source_exposure_length": str(row.source_exposure_length),
                "allocated_transition_weight": str(row.allocated_transition_weight),
                "time_stratified_rate": str(row.time_stratified_rate),
                "global_rate": str(row.global_rate),
            }
            for row in report.matrix_rows
        ],
    )


def write_time_stratified_branch_table(
    path: Path,
    report: TimeStratifiedTransitionReport,
) -> Path:
    """Write one branch-overlap ledger across time bins."""
    return write_taxon_rows(
        path,
        columns=[
            "time_bin_label",
            "start_depth",
            "end_depth",
            "parent_node",
            "child_node",
            "parent_depth",
            "child_depth",
            "source_region",
            "target_region",
            "changed",
            "overlap_length",
            "allocated_transition_weight",
            "support",
            "strongly_supported",
        ],
        rows=[
            {
                "time_bin_label": row.time_bin_label,
                "start_depth": str(row.start_depth),
                "end_depth": str(row.end_depth),
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "parent_depth": str(row.parent_depth),
                "child_depth": str(row.child_depth),
                "source_region": row.source_region,
                "target_region": row.target_region,
                "changed": str(row.changed).lower(),
                "overlap_length": str(row.overlap_length),
                "allocated_transition_weight": str(row.allocated_transition_weight),
                "support": str(row.support),
                "strongly_supported": str(row.strongly_supported).lower(),
            }
            for row in report.branch_rows
        ],
    )


def write_time_stratified_exclusion_table(
    path: Path,
    report: TimeStratifiedTransitionReport,
) -> Path:
    """Write one excluded-taxa ledger for time-stratified geography."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "raw_region",
            "normalized_region",
            "reason",
            "note",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "raw_region": row.raw_region,
                "normalized_region": row.normalized_region or "",
                "reason": row.reason,
                "note": row.note,
            }
            for row in report.exclusion_rows
        ],
    )


def _normalize_time_bins(time_bins: list[TimeBinDefinition]) -> list[TimeBinDefinition]:
    if not time_bins:
        raise ValueError(
            "time-stratified geographic transitions require at least one time bin"
        )
    normalized = sorted(
        time_bins, key=lambda row: (row.start_depth, row.end_depth, row.label)
    )
    seen_labels: set[str] = set()
    previous_end = 0.0
    for index, row in enumerate(normalized):
        if not row.label.strip():
            raise ValueError("time bin labels must be non-empty")
        if row.label in seen_labels:
            raise ValueError(f"duplicate time bin label: {row.label}")
        seen_labels.add(row.label)
        if row.start_depth < 0.0:
            raise ValueError(f"time bin {row.label} starts before the root")
        if row.end_depth <= row.start_depth:
            raise ValueError(f"time bin {row.label} must end after it starts")
        if index and row.start_depth < previous_end:
            raise ValueError("time bins must not overlap")
        previous_end = row.end_depth
    return normalized


def _validate_time_stratified_tree(tree, tree_depth: float) -> None:
    branch_lengths = tree.branch_lengths()
    if not branch_lengths or any(length is None for length in branch_lengths):
        raise AncestralReconstructionError(
            "time-stratified geographic transitions require complete branch lengths"
        )
    if any((length or 0.0) <= 0.0 for length in branch_lengths):
        raise AncestralReconstructionError(
            "time-stratified geographic transitions require strictly positive branch lengths"
        )
    if tree_depth <= 0.0:
        raise AncestralReconstructionError(
            "time-stratified geographic transitions require positive tree depth"
        )


def _validate_time_bins(time_bins: list[TimeBinDefinition], tree_depth: float) -> None:
    for row in time_bins:
        if row.end_depth > tree_depth + 1e-12:
            raise ValueError(
                f"time bin {row.label} extends beyond tree depth {tree_depth}"
            )


def _build_time_bin_coverage_warning(
    time_bins: list[TimeBinDefinition],
    tree_depth: float,
) -> str | None:
    tolerance = 1e-12
    uncovered_ranges: list[tuple[float, float]] = []
    cursor = 0.0
    for row in time_bins:
        if row.start_depth > cursor + tolerance:
            uncovered_ranges.append((cursor, row.start_depth))
        cursor = max(cursor, row.end_depth)
    if cursor < tree_depth - tolerance:
        uncovered_ranges.append((cursor, tree_depth))
    if not uncovered_ranges:
        return None
    formatted_ranges = ", ".join(
        f"[{stable_value(start)}, {stable_value(end)})"
        for start, end in uncovered_ranges
    )
    return (
        "time bins do not cover the full tree depth; branch segments outside the "
        f"requested intervals are excluded from interval-specific rates: "
        f"{formatted_ranges}"
    )


def _node_depths(tree) -> dict[str, float]:
    depths: dict[str, float] = {node_signature(tree.root): 0.0}

    def visit(node, depth: float) -> None:
        current_signature = node_signature(node)
        depths[current_signature] = stable_value(depth)
        for child in node.children:
            child_depth = depth + float(child.branch_length or 0.0)
            visit(child, child_depth)

    visit(tree.root, 0.0)
    return depths


def _tree_depth(tree) -> float:
    lengths = [length for length in tree.root_to_tip_lengths() if length is not None]
    if not lengths:
        return 0.0
    return stable_value(max(float(length) for length in lengths))


def _build_branch_rows(
    base_report: GeographicStateModelReport,
    time_bins: list[TimeBinDefinition],
    depth_by_node: dict[str, float],
) -> list[TimeStratifiedBranchRow]:
    rows: list[TimeStratifiedBranchRow] = []
    for event in base_report.transition_event_rows:
        parent_depth = depth_by_node[event.parent_node]
        child_depth = depth_by_node[event.child_node]
        branch_length = child_depth - parent_depth
        for time_bin in time_bins:
            overlap = _overlap_length(
                parent_depth,
                child_depth,
                time_bin.start_depth,
                time_bin.end_depth,
            )
            if overlap <= 0.0:
                continue
            allocated_weight = (
                stable_value(overlap / branch_length) if event.changed else 0.0
            )
            rows.append(
                TimeStratifiedBranchRow(
                    time_bin_label=time_bin.label,
                    start_depth=time_bin.start_depth,
                    end_depth=time_bin.end_depth,
                    parent_node=event.parent_node,
                    child_node=event.child_node,
                    parent_depth=parent_depth,
                    child_depth=child_depth,
                    source_region=event.source_region,
                    target_region=event.target_region,
                    changed=event.changed,
                    overlap_length=stable_value(overlap),
                    allocated_transition_weight=allocated_weight,
                    support=stable_value(event.support),
                    strongly_supported=event.strongly_supported,
                )
            )
    return rows


def _build_matrix_rows(
    base_report: GeographicStateModelReport,
    time_bins: list[TimeBinDefinition],
    branch_rows: list[TimeStratifiedBranchRow],
) -> list[TimeStratifiedTransitionMatrixRow]:
    global_rate_by_pair = {
        (row.source_region, row.target_region): row.rate
        for row in base_report.transition_rate_rows
    }
    source_exposure: dict[tuple[str, str], float] = {}
    allocated_change: dict[tuple[str, str, str], float] = {}
    for row in branch_rows:
        source_key = (row.time_bin_label, row.source_region)
        source_exposure[source_key] = (
            source_exposure.get(source_key, 0.0) + row.overlap_length
        )
        change_key = (row.time_bin_label, row.source_region, row.target_region)
        allocated_change[change_key] = (
            allocated_change.get(change_key, 0.0) + row.allocated_transition_weight
        )
    rows: list[TimeStratifiedTransitionMatrixRow] = []
    for time_bin in time_bins:
        for source_region, target_region in sorted(global_rate_by_pair):
            exposure = stable_value(
                source_exposure.get((time_bin.label, source_region), 0.0)
            )
            change_weight = stable_value(
                allocated_change.get(
                    (time_bin.label, source_region, target_region), 0.0
                )
            )
            rate = 0.0 if exposure == 0.0 else stable_value(change_weight / exposure)
            rows.append(
                TimeStratifiedTransitionMatrixRow(
                    time_bin_label=time_bin.label,
                    start_depth=time_bin.start_depth,
                    end_depth=time_bin.end_depth,
                    source_region=source_region,
                    target_region=target_region,
                    source_exposure_length=exposure,
                    allocated_transition_weight=change_weight,
                    time_stratified_rate=rate,
                    global_rate=stable_value(
                        global_rate_by_pair[(source_region, target_region)]
                    ),
                )
            )
    return rows


def _overlap_length(
    left_start: float,
    left_end: float,
    right_start: float,
    right_end: float,
) -> float:
    return max(0.0, min(left_end, right_end) - max(left_start, right_start))
