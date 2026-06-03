from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path
from statistics import mean

from bijux_phylogenetics.datasets.study_inputs import (
    align_tree_and_trait_table,
    validate_traits_table,
    write_taxon_rows,
)
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.node_identity import build_ape_internal_node_map
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from ..common import node_signature, tip_root_depths


@dataclass(slots=True)
class DisparityTaxonExclusion:
    """One taxon excluded before continuous disparity analysis."""

    taxon: str
    reason: str


@dataclass(slots=True)
class CladeDisparityRow:
    """One internal clade with geiger-style continuous disparity evidence."""

    ape_node_id: int
    node_kind: str
    node_label: str | None
    node: str
    descendant_taxa: list[str]
    descendant_taxon_count: int
    branch_length_depth: float
    branching_time: float
    relative_branching_time: float
    stem_branching_time: float
    relative_stem_branching_time: float
    disparity: float


@dataclass(slots=True)
class ContinuousCladeDisparityReport:
    """Clade-wise continuous disparity surface for one rooted ultrametric tree."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait_columns: list[str]
    tree_taxon_count: int
    analyzed_taxa: list[str]
    excluded_taxa: list[DisparityTaxonExclusion]
    distance_metric: str
    method_formula: str
    root_age: float
    root_disparity: float
    minimum_clade_disparity: float
    maximum_clade_disparity: float
    clade_rows: list[CladeDisparityRow]
    warnings: list[str]
    assumptions: list[str]

    @property
    def analyzed_taxon_count(self) -> int:
        return len(self.analyzed_taxa)


@dataclass(slots=True)
class DisparityThroughTimeCurveRow:
    """One geiger-style disparity-through-time curve point."""

    point_index: int
    relative_time: float
    raw_disparity: float
    relative_disparity: float
    contributing_clade_count: int


@dataclass(slots=True)
class DisparityThroughTimeBinRow:
    """One optional equal-width time bin over the DTT curve."""

    bin_index: int
    start_time: float
    end_time: float
    midpoint_time: float
    point_count: int
    mean_relative_disparity: float | None
    minimum_relative_disparity: float | None
    maximum_relative_disparity: float | None


@dataclass(slots=True)
class DisparityThroughTimeReport:
    """Reviewer-facing disparity-through-time report aligned to geiger DTT."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait_columns: list[str]
    tree_taxon_count: int
    analyzed_taxa: list[str]
    excluded_taxa: list[DisparityTaxonExclusion]
    distance_metric: str
    root_age: float
    root_disparity: float
    relative_scaling_applied: bool
    clade_rows: list[CladeDisparityRow]
    curve_rows: list[DisparityThroughTimeCurveRow]
    time_bin_rows: list[DisparityThroughTimeBinRow]
    warnings: list[str]
    assumptions: list[str]

    @property
    def analyzed_taxon_count(self) -> int:
        return len(self.analyzed_taxa)


def summarize_continuous_clade_disparity(
    tree_path: Path,
    traits_path: Path,
    *,
    trait_columns: list[str],
    taxon_column: str | None = None,
) -> ContinuousCladeDisparityReport:
    """Summarize geiger-style continuous disparity for every internal clade."""
    selected_traits = _normalize_trait_columns(trait_columns)
    dataset = _load_continuous_disparity_dataset(
        tree_path,
        traits_path,
        trait_columns=selected_traits,
        taxon_column=taxon_column,
    )
    clade_rows = _build_clade_disparity_rows(
        dataset.tree,
        dataset.root_age,
        dataset.values_by_taxon,
        analyzed_taxa=dataset.analyzed_taxa,
    )
    return ContinuousCladeDisparityReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait_columns=selected_traits,
        tree_taxon_count=dataset.tree_taxon_count,
        analyzed_taxa=dataset.analyzed_taxa,
        excluded_taxa=dataset.excluded_taxa,
        distance_metric="avg-squared-euclidean",
        method_formula=(
            "mean over all tip pairs of the squared Euclidean distance across the requested continuous trait matrix"
        ),
        root_age=dataset.root_age,
        root_disparity=clade_rows[0].disparity,
        minimum_clade_disparity=min(row.disparity for row in clade_rows),
        maximum_clade_disparity=max(row.disparity for row in clade_rows),
        clade_rows=clade_rows,
        warnings=dataset.warnings,
        assumptions=[
            "Clade disparity matches geiger avg.sq disparity, defined as the mean squared Euclidean distance across all trait-vector pairs inside each internal clade.",
            "Continuous disparity is computed after explicit tree-versus-trait alignment and pruning of overlapping taxa with missing or non-numeric values in any requested trait column.",
            "Internal clades are reported in ape-style preorder with the root first so downstream disparity-through-time summaries can match geiger node ordering.",
        ],
    )


def write_continuous_clade_disparity_summary_table(
    path: Path,
    report: ContinuousCladeDisparityReport,
) -> Path:
    """Write one one-row summary ledger for continuous clade disparity."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon_column",
            "trait_columns",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "clade_count",
            "distance_metric",
            "method_formula",
            "root_age",
            "root_disparity",
            "minimum_clade_disparity",
            "maximum_clade_disparity",
            "warning_count",
        ],
        rows=[
            {
                "taxon_column": report.taxon_column,
                "trait_columns": "|".join(report.trait_columns),
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": report.analyzed_taxon_count,
                "excluded_taxon_count": len(report.excluded_taxa),
                "clade_count": len(report.clade_rows),
                "distance_metric": report.distance_metric,
                "method_formula": report.method_formula,
                "root_age": report.root_age,
                "root_disparity": report.root_disparity,
                "minimum_clade_disparity": report.minimum_clade_disparity,
                "maximum_clade_disparity": report.maximum_clade_disparity,
                "warning_count": len(report.warnings),
            }
        ],
    )


def summarize_disparity_through_time(
    tree_path: Path,
    traits_path: Path,
    *,
    trait_columns: list[str],
    taxon_column: str | None = None,
    time_bin_count: int | None = None,
) -> DisparityThroughTimeReport:
    """Summarize geiger-style disparity through time for continuous traits."""
    if time_bin_count is not None and time_bin_count < 1:
        raise ComparativeMethodError(
            "disparity-through-time time-bin count must be at least one"
        )
    clade_report = summarize_continuous_clade_disparity(
        tree_path,
        traits_path,
        trait_columns=trait_columns,
        taxon_column=taxon_column,
    )
    root_disparity = clade_report.root_disparity
    normalized_branching_times = sorted(
        (row.relative_branching_time for row in clade_report.clade_rows),
        reverse=True,
    )
    raw_curve = [root_disparity]
    contributing_counts = [1]
    for threshold in normalized_branching_times[:-1]:
        contributing = [
            row.disparity
            for row in clade_report.clade_rows
            if row.relative_stem_branching_time >= threshold
            and row.relative_branching_time < threshold
        ]
        raw_curve.append(mean(contributing) if contributing else 0.0)
        contributing_counts.append(len(contributing))
    raw_curve.append(0.0)
    contributing_counts.append(0)
    curve_times = [0.0] + [1.0 - value for value in normalized_branching_times]
    relative_scaling_applied = root_disparity > 0.0
    curve_rows = [
        DisparityThroughTimeCurveRow(
            point_index=index + 1,
            relative_time=float(format(relative_time, ".15g")),
            raw_disparity=float(format(raw_value, ".15g")),
            relative_disparity=float(
                format(
                    (raw_value / root_disparity)
                    if relative_scaling_applied
                    else raw_value,
                    ".15g",
                )
            ),
            contributing_clade_count=contributing_counts[index],
        )
        for index, (relative_time, raw_value) in enumerate(
            zip(curve_times, raw_curve, strict=True)
        )
    ]
    time_bin_rows = _build_time_bin_rows(curve_rows, time_bin_count=time_bin_count)
    warnings = list(clade_report.warnings)
    if not relative_scaling_applied:
        warnings.append(
            "root clade disparity is zero, so the disparity-through-time curve remains on the raw disparity scale instead of relative disparity"
        )
    return DisparityThroughTimeReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=clade_report.taxon_column,
        trait_columns=list(clade_report.trait_columns),
        tree_taxon_count=clade_report.tree_taxon_count,
        analyzed_taxa=list(clade_report.analyzed_taxa),
        excluded_taxa=list(clade_report.excluded_taxa),
        distance_metric=clade_report.distance_metric,
        root_age=clade_report.root_age,
        root_disparity=float(format(root_disparity, ".15g")),
        relative_scaling_applied=relative_scaling_applied,
        clade_rows=clade_report.clade_rows,
        curve_rows=curve_rows,
        time_bin_rows=time_bin_rows,
        warnings=warnings,
        assumptions=[
            "Disparity-through-time follows geiger DTT by averaging clade disparities across internal clades whose stems cross each successive branching-time threshold.",
            "Reported times run from root to present on a relative 0-to-1 axis, matching the geiger DTT output ordering with an explicit leading root point.",
            *clade_report.assumptions,
        ],
    )


def write_continuous_clade_disparity_table(
    path: Path,
    report: ContinuousCladeDisparityReport,
) -> Path:
    """Write one long-form clade disparity ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "ape_node_id",
            "node_kind",
            "node_label",
            "node",
            "descendant_taxa",
            "descendant_taxon_count",
            "branch_length_depth",
            "branching_time",
            "relative_branching_time",
            "stem_branching_time",
            "relative_stem_branching_time",
            "disparity",
            "trait_columns",
            "distance_metric",
        ],
        rows=[
            {
                "ape_node_id": row.ape_node_id,
                "node_kind": row.node_kind,
                "node_label": row.node_label or "",
                "node": row.node,
                "descendant_taxa": "|".join(row.descendant_taxa),
                "descendant_taxon_count": row.descendant_taxon_count,
                "branch_length_depth": row.branch_length_depth,
                "branching_time": row.branching_time,
                "relative_branching_time": row.relative_branching_time,
                "stem_branching_time": row.stem_branching_time,
                "relative_stem_branching_time": row.relative_stem_branching_time,
                "disparity": row.disparity,
                "trait_columns": "|".join(report.trait_columns),
                "distance_metric": report.distance_metric,
            }
            for row in report.clade_rows
        ],
    )


def write_disparity_through_time_summary_table(
    path: Path,
    report: DisparityThroughTimeReport,
) -> Path:
    """Write one one-row DTT summary ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon_column",
            "trait_columns",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "clade_count",
            "curve_point_count",
            "time_bin_count",
            "distance_metric",
            "root_age",
            "root_disparity",
            "relative_scaling_applied",
            "warning_count",
        ],
        rows=[
            {
                "taxon_column": report.taxon_column,
                "trait_columns": "|".join(report.trait_columns),
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": report.analyzed_taxon_count,
                "excluded_taxon_count": len(report.excluded_taxa),
                "clade_count": len(report.clade_rows),
                "curve_point_count": len(report.curve_rows),
                "time_bin_count": len(report.time_bin_rows),
                "distance_metric": report.distance_metric,
                "root_age": report.root_age,
                "root_disparity": report.root_disparity,
                "relative_scaling_applied": report.relative_scaling_applied,
                "warning_count": len(report.warnings),
            }
        ],
    )


def write_disparity_through_time_curve_table(
    path: Path,
    report: DisparityThroughTimeReport,
) -> Path:
    """Write one long-form DTT curve ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "point_index",
            "relative_time",
            "raw_disparity",
            "relative_disparity",
            "contributing_clade_count",
            "trait_columns",
            "distance_metric",
        ],
        rows=[
            {
                "point_index": row.point_index,
                "relative_time": row.relative_time,
                "raw_disparity": row.raw_disparity,
                "relative_disparity": row.relative_disparity,
                "contributing_clade_count": row.contributing_clade_count,
                "trait_columns": "|".join(report.trait_columns),
                "distance_metric": report.distance_metric,
            }
            for row in report.curve_rows
        ],
    )


def write_disparity_through_time_bin_table(
    path: Path,
    report: DisparityThroughTimeReport,
) -> Path:
    """Write one optional equal-width time-bin ledger over the DTT curve."""
    return write_taxon_rows(
        path,
        columns=[
            "bin_index",
            "start_time",
            "end_time",
            "midpoint_time",
            "point_count",
            "mean_relative_disparity",
            "minimum_relative_disparity",
            "maximum_relative_disparity",
            "trait_columns",
        ],
        rows=[
            {
                "bin_index": row.bin_index,
                "start_time": row.start_time,
                "end_time": row.end_time,
                "midpoint_time": row.midpoint_time,
                "point_count": row.point_count,
                "mean_relative_disparity": (
                    ""
                    if row.mean_relative_disparity is None
                    else row.mean_relative_disparity
                ),
                "minimum_relative_disparity": (
                    ""
                    if row.minimum_relative_disparity is None
                    else row.minimum_relative_disparity
                ),
                "maximum_relative_disparity": (
                    ""
                    if row.maximum_relative_disparity is None
                    else row.maximum_relative_disparity
                ),
                "trait_columns": "|".join(report.trait_columns),
            }
            for row in report.time_bin_rows
        ],
    )


def write_disparity_through_time_exclusion_table(
    path: Path,
    report: DisparityThroughTimeReport,
) -> Path:
    """Write one excluded-taxa ledger for DTT."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason", "trait_columns"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
                "trait_columns": "|".join(report.trait_columns),
            }
            for row in report.excluded_taxa
        ],
    )


def render_disparity_through_time_svg(
    path: Path,
    report: DisparityThroughTimeReport,
) -> int:
    """Render one reviewer-facing SVG DTT curve artifact."""
    width = 960
    height = 360
    left = 92
    right = 36
    top = 34
    bottom = 58
    plot_width = width - left - right
    plot_height = height - top - bottom
    maximum_disparity = max(
        max((row.relative_disparity for row in report.curve_rows), default=0.0),
        1e-9,
    )
    x_ticks = [0.0, 0.5, 1.0]
    y_ticks = [0.0, 0.25, 0.5, 0.75, 1.0]

    def x_position(relative_time: float) -> float:
        return left + relative_time * plot_width

    def y_position(relative_disparity: float) -> float:
        return top + (1.0 - (relative_disparity / maximum_disparity)) * plot_height

    segments = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="disparity through time curve">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#f8fafc" />',
        f'<text x="{left}" y="18" font-size="18" font-family="Avenir Next, Segoe UI, sans-serif" fill="#0f172a">Disparity-through-time curve</text>',
        f'<text x="{left}" y="{height - 18}" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#475569">relative time from root to present</text>',
        f'<text x="18" y="{top - 10}" font-size="12" font-family="Avenir Next, Segoe UI, sans-serif" fill="#475569">relative disparity</text>',
        f'<text x="{left}" y="{height - 34}" font-size="11" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#64748b">{escape(", ".join(report.trait_columns))}</text>',
        f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#cbd5e1" stroke-width="1.5" />',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#cbd5e1" stroke-width="1.5" />',
    ]
    for tick in x_ticks:
        x = x_position(tick)
        segments.extend(
            [
                f'<line x1="{x}" y1="{top}" x2="{x}" y2="{top + plot_height}" stroke="#e2e8f0" stroke-width="1" />',
                f'<text x="{x}" y="{top + plot_height + 22}" text-anchor="middle" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#475569">{_format_number(tick)}</text>',
            ]
        )
    for tick in y_ticks:
        y = y_position(tick * maximum_disparity)
        segments.extend(
            [
                f'<line x1="{left}" y1="{y}" x2="{left + plot_width}" y2="{y}" stroke="#e2e8f0" stroke-width="1" />',
                f'<text x="{left - 14}" y="{y + 4}" text-anchor="end" font-size="12" font-family="Iowan Old Style, Palatino Linotype, serif" fill="#475569">{_format_number(tick * maximum_disparity if maximum_disparity != 1.0 else tick)}</text>',
            ]
        )
    polyline = " ".join(
        f"{x_position(row.relative_time):.2f},{y_position(row.relative_disparity):.2f}"
        for row in report.curve_rows
    )
    segments.append(
        f'<polyline points="{polyline}" fill="none" stroke="#0f766e" stroke-width="4" stroke-linejoin="round" stroke-linecap="round" />'
    )
    for row in report.curve_rows:
        x = x_position(row.relative_time)
        y = y_position(row.relative_disparity)
        segments.append(
            f'<circle cx="{x}" cy="{y}" r="4.5" fill="#0f766e" stroke="#ffffff" stroke-width="1.5" />'
        )
    segments.append("</svg>")
    path.write_text("\n".join(segments) + "\n", encoding="utf-8")
    return len(report.curve_rows)


@dataclass(slots=True)
class _ContinuousDisparityDataset:
    tree: PhyloTree
    taxon_column: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    values_by_taxon: dict[str, list[float]]
    excluded_taxa: list[DisparityTaxonExclusion]
    root_age: float
    warnings: list[str]


def _normalize_trait_columns(trait_columns: list[str]) -> list[str]:
    normalized: list[str] = []
    for column in trait_columns:
        if column not in normalized:
            normalized.append(column)
    if not normalized:
        raise ComparativeMethodError(
            "disparity-through-time analysis requires at least one trait column"
        )
    return normalized


def _load_continuous_disparity_dataset(
    tree_path: Path,
    traits_path: Path,
    *,
    trait_columns: list[str],
    taxon_column: str | None,
) -> _ContinuousDisparityDataset:
    validation = validate_traits_table(traits_path, taxon_column=taxon_column)
    by_name = {column.name: column for column in validation.trait_columns}
    for trait in trait_columns:
        summary = by_name.get(trait)
        if summary is None:
            raise ComparativeMethodError(
                f"trait table does not contain column '{trait}'"
            )
        if summary.kind != "numeric":
            raise ComparativeMethodError(
                f"trait column '{trait}' must be numeric for disparity-through-time analysis"
            )
    alignment = align_tree_and_trait_table(
        tree_path,
        traits_path,
        taxon_column=taxon_column,
        required_trait_columns=trait_columns,
        drop_missing_for_columns=trait_columns,
    )
    rows_by_taxon = {row[alignment.taxon_column]: row for row in alignment.rows}
    analyzed_taxa: list[str] = []
    values_by_taxon: dict[str, list[float]] = {}
    excluded_taxa: list[DisparityTaxonExclusion] = []
    non_numeric_taxa: list[str] = []
    for taxon in alignment.report.aligned_taxa:
        row = rows_by_taxon[taxon]
        try:
            values = [float(row[trait]) for trait in trait_columns]
        except ValueError:
            non_numeric_taxa.append(taxon)
            continue
        analyzed_taxa.append(taxon)
        values_by_taxon[taxon] = values
    excluded_taxa.extend(
        DisparityTaxonExclusion(taxon=taxon, reason="missing_from_trait_table")
        for taxon in alignment.report.dropped_tree_taxa
    )
    excluded_taxa.extend(
        DisparityTaxonExclusion(taxon=taxon, reason="absent_from_tree")
        for taxon in alignment.report.dropped_trait_taxa
    )
    excluded_taxa.extend(
        DisparityTaxonExclusion(taxon=taxon, reason="missing_trait_value")
        for taxon in alignment.report.dropped_missing_value_taxa
    )
    excluded_taxa.extend(
        DisparityTaxonExclusion(taxon=taxon, reason="non_numeric_trait_value")
        for taxon in non_numeric_taxa
    )
    if len(analyzed_taxa) < 3:
        raise ComparativeMethodError(
            "disparity-through-time analysis requires at least three taxa with numeric trait values"
        )
    if non_numeric_taxa:
        pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, analyzed_taxa)
    else:
        pruned_tree = alignment.tree
    _validate_dtt_tree(pruned_tree)
    root_age = _root_age(pruned_tree)
    warnings: list[str] = []
    exclusion_reasons = {row.reason for row in excluded_taxa}
    if "missing_from_trait_table" in exclusion_reasons:
        warnings.append(
            "trait table is missing one or more tree taxa and those taxa were pruned"
        )
    if "absent_from_tree" in exclusion_reasons:
        warnings.append("trait table contains taxa absent from the tree")
    if "missing_trait_value" in exclusion_reasons:
        warnings.append(
            "one or more overlapping taxa have missing trait values and were pruned"
        )
    if "non_numeric_trait_value" in exclusion_reasons:
        warnings.append(
            "one or more overlapping taxa have non-numeric trait values and were pruned"
        )
    return _ContinuousDisparityDataset(
        tree=pruned_tree,
        taxon_column=alignment.taxon_column,
        tree_taxon_count=alignment.report.original_tree_taxa,
        analyzed_taxa=analyzed_taxa,
        values_by_taxon=values_by_taxon,
        excluded_taxa=excluded_taxa,
        root_age=root_age,
        warnings=warnings,
    )


def _validate_dtt_tree(tree: PhyloTree) -> None:
    if len(tree.root.children) != 2:
        raise ComparativeMethodError(
            "disparity-through-time analysis requires a rooted tree"
        )
    minimum_branch_length = min(
        float(node.branch_length or 0.0)
        for node in tree.iter_nodes()
        if node is not tree.root
    )
    if any(
        node.branch_length is None
        for node in tree.iter_nodes()
        if node is not tree.root
    ):
        raise ComparativeMethodError(
            "disparity-through-time analysis requires complete branch lengths"
        )
    if minimum_branch_length < 0.0:
        raise ComparativeMethodError(
            "disparity-through-time analysis requires non-negative branch lengths"
        )
    tip_depths = tip_root_depths(tree, tree.tip_names)
    maximum_tip_depth = max(tip_depths.values())
    minimum_tip_depth = min(tip_depths.values())
    if maximum_tip_depth - minimum_tip_depth > 1e-12:
        raise ComparativeMethodError(
            "disparity-through-time analysis requires an ultrametric tree"
        )
    if maximum_tip_depth <= 0.0:
        raise ComparativeMethodError(
            "disparity-through-time analysis requires a positive root age"
        )


def _root_age(tree: PhyloTree) -> float:
    return max(tip_root_depths(tree, tree.tip_names).values())


def _build_clade_disparity_rows(
    tree: PhyloTree,
    root_age: float,
    values_by_taxon: dict[str, list[float]],
    *,
    analyzed_taxa: list[str],
) -> list[CladeDisparityRow]:
    depth_lookup = _node_depth_lookup(tree)
    rows: list[CladeDisparityRow] = []
    for ape_node_id, node in build_ape_internal_node_map(tree).items():
        descendant_taxa = [
            taxon for taxon in analyzed_taxa if taxon in set(node.descendant_taxa)
        ]
        matrix = [values_by_taxon[taxon] for taxon in descendant_taxa]
        branch_length_depth = depth_lookup[node.node_id or ""]
        branching_time = root_age - branch_length_depth
        parent_depth = (
            branch_length_depth
            if node is tree.root
            else depth_lookup[node.parent.node_id or ""]
        )
        stem_branching_time = root_age - parent_depth
        rows.append(
            CladeDisparityRow(
                ape_node_id=ape_node_id,
                node_kind="root" if node is tree.root else "internal",
                node_label=node.name,
                node=node_signature(node),
                descendant_taxa=descendant_taxa,
                descendant_taxon_count=len(descendant_taxa),
                branch_length_depth=float(format(branch_length_depth, ".15g")),
                branching_time=float(format(branching_time, ".15g")),
                relative_branching_time=float(
                    format(branching_time / root_age, ".15g")
                ),
                stem_branching_time=float(format(stem_branching_time, ".15g")),
                relative_stem_branching_time=float(
                    format(stem_branching_time / root_age, ".15g")
                ),
                disparity=float(
                    format(_average_squared_euclidean_distance(matrix), ".15g")
                ),
            )
        )
    return rows


def _average_squared_euclidean_distance(matrix: list[list[float]]) -> float:
    if len(matrix) < 2:
        return 0.0
    distances: list[float] = []
    for left_index in range(len(matrix) - 1):
        left = matrix[left_index]
        for right in matrix[left_index + 1 :]:
            distances.append(
                sum(
                    (left_value - right_value) ** 2
                    for left_value, right_value in zip(left, right, strict=True)
                )
            )
    return mean(distances) if distances else 0.0


def _build_time_bin_rows(
    curve_rows: list[DisparityThroughTimeCurveRow],
    *,
    time_bin_count: int | None,
) -> list[DisparityThroughTimeBinRow]:
    if time_bin_count is None:
        return []
    bin_width = 1.0 / time_bin_count
    rows: list[DisparityThroughTimeBinRow] = []
    for index in range(time_bin_count):
        start = index * bin_width
        end = 1.0 if index == time_bin_count - 1 else (index + 1) * bin_width
        values = [
            row.relative_disparity
            for row in curve_rows
            if (start <= row.relative_time < end)
            or (index == time_bin_count - 1 and row.relative_time == end)
        ]
        rows.append(
            DisparityThroughTimeBinRow(
                bin_index=index + 1,
                start_time=float(format(start, ".15g")),
                end_time=float(format(end, ".15g")),
                midpoint_time=float(format((start + end) / 2.0, ".15g")),
                point_count=len(values),
                mean_relative_disparity=(
                    None if not values else float(format(mean(values), ".15g"))
                ),
                minimum_relative_disparity=(
                    None if not values else float(format(min(values), ".15g"))
                ),
                maximum_relative_disparity=(
                    None if not values else float(format(max(values), ".15g"))
                ),
            )
        )
    return rows


def _node_depth_lookup(tree: PhyloTree) -> dict[str, float]:
    depths: dict[str, float] = {tree.root.node_id or "": 0.0}

    def visit(node: TreeNode) -> None:
        node_id = node.node_id or ""
        current_depth = depths[node_id]
        for child in node.children:
            depths[child.node_id or ""] = current_depth + float(
                child.branch_length or 0.0
            )
            visit(child)

    visit(tree.root)
    return depths


def _format_number(value: float) -> str:
    if abs(value) < 1e-12:
        return "0"
    return format(value, ".3g")
