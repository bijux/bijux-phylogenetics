from __future__ import annotations

from html import escape
from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylogeography.continuous_coordinate_review.artifact_outputs import (
    write_coordinate_estimate_table,
    write_coordinate_movement_branch_table,
    write_coordinate_movement_exclusion_table,
    write_coordinate_movement_outlier_table,
    write_coordinate_movement_summary_table,
)
from bijux_phylogenetics.phylogeography.continuous_coordinate_review.contracts import (
    CoordinateEstimateRow,
    CoordinateMovementBranchRow,
    CoordinateMovementExclusionRow,
    CoordinateMovementOutlierRow,
    CoordinateMovementSummary,
    CoordinateMovementVisualization,
    PhylogeographicCoordinateReport,
)
from bijux_phylogenetics.phylogeography.continuous_coordinate_review.inputs import (
    prepare_coordinate_dataset as _prepare_coordinate_dataset,
    write_filtered_coordinate_table as _write_filtered_coordinate_table,
)
from bijux_phylogenetics.phylogeography.continuous_coordinate_review.movement_analysis import (
    build_branch_rows as _build_branch_rows,
    build_estimate_rows as _build_estimate_rows,
    build_summary as _build_summary,
    review_branch_outliers as _review_branch_outliers,
    stable_float as _stable_float,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

__all__ = [
    "CoordinateEstimateRow",
    "CoordinateMovementBranchRow",
    "CoordinateMovementExclusionRow",
    "CoordinateMovementOutlierRow",
    "CoordinateMovementSummary",
    "CoordinateMovementVisualization",
    "PhylogeographicCoordinateReport",
    "render_coordinate_movement_visualization",
    "summarize_continuous_phylogeography",
    "write_coordinate_estimate_table",
    "write_coordinate_movement_branch_table",
    "write_coordinate_movement_exclusion_table",
    "write_coordinate_movement_outlier_table",
    "write_coordinate_movement_summary_table",
]


def summarize_continuous_phylogeography(
    tree_path: Path,
    table_path: Path,
    *,
    latitude_column: str,
    longitude_column: str,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
) -> PhylogeographicCoordinateReport:
    """Reconstruct ancestral coordinates and review branchwise movement."""
    prepared = _prepare_coordinate_dataset(
        tree_path,
        table_path,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
        taxon_column=taxon_column,
    )
    if len(prepared.included_taxa) < 2:
        raise AncestralReconstructionError(
            "continuous phylogeography requires at least two taxa with usable coordinates"
        )

    with tempfile.TemporaryDirectory(prefix="bijux-phylogeography-") as temp_dir:
        temp_root = Path(temp_dir)
        pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, prepared.included_taxa)
        pruned_tree_path = temp_root / "analysis-tree.nwk"
        pruned_tree_path.write_text(f"{dumps_newick(pruned_tree)}\n", encoding="utf-8")
        filtered_table_path = temp_root / "coordinates.tsv"
        _write_filtered_coordinate_table(filtered_table_path, prepared)
        latitude_report = reconstruct_continuous_ancestral_states(
            pruned_tree_path,
            filtered_table_path,
            trait=latitude_column,
            taxon_column=prepared.taxon_column,
            model=model,
            alpha=alpha,
        )
        longitude_report = reconstruct_continuous_ancestral_states(
            pruned_tree_path,
            filtered_table_path,
            trait=longitude_column,
            taxon_column=prepared.taxon_column,
            model=model,
            alpha=alpha,
        )

    warnings = list(
        dict.fromkeys([*latitude_report.warnings, *longitude_report.warnings])
    )
    estimate_rows = _build_estimate_rows(latitude_report, longitude_report)
    branch_rows = _build_branch_rows(latitude_report, longitude_report)
    outlier_rows, reviewed_branch_rows = _review_branch_outliers(branch_rows)
    summary = _build_summary(
        report=latitude_report,
        prepared=prepared,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
        estimate_rows=estimate_rows,
        branch_rows=reviewed_branch_rows,
        outlier_rows=outlier_rows,
        warnings=warnings,
    )
    return PhylogeographicCoordinateReport(
        tree_path=tree_path,
        table_path=table_path,
        taxon_column=prepared.taxon_column,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
        model=model,
        alpha=alpha,
        summary=summary,
        estimate_rows=estimate_rows,
        branch_rows=reviewed_branch_rows,
        outlier_rows=outlier_rows,
        exclusion_rows=prepared.exclusion_rows,
        warnings=warnings,
    )


def render_coordinate_movement_visualization(
    report: PhylogeographicCoordinateReport,
    *,
    out_path: Path,
    width: int = 960,
    height: int = 640,
) -> CoordinateMovementVisualization:
    """Render a coordinate-space movement visualization as SVG or HTML."""
    output_format = out_path.suffix.lower().lstrip(".")
    if output_format not in {"svg", "html"}:
        raise ValueError(
            "coordinate movement visualization output must end in .svg or .html"
        )
    svg_markup = _build_coordinate_svg(report, width=width, height=height)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if output_format == "svg":
        out_path.write_text(svg_markup, encoding="utf-8")
    else:
        out_path.write_text(
            _build_coordinate_html(report, svg_markup), encoding="utf-8"
        )
    return CoordinateMovementVisualization(
        output_path=out_path,
        format=output_format,
        width=width,
        height=height,
        highlighted_branch_count=sum(
            row.outlier_jump or row.impossible_jump for row in report.branch_rows
        ),
    )


def _build_coordinate_svg(
    report: PhylogeographicCoordinateReport,
    *,
    width: int,
    height: int,
) -> str:
    rows = report.estimate_rows
    latitudes = [row.latitude for row in rows]
    longitudes = [row.longitude for row in rows]
    min_latitude = min(latitudes)
    max_latitude = max(latitudes)
    min_longitude = min(longitudes)
    max_longitude = max(longitudes)
    padding = 48.0

    def project(longitude: float, latitude: float) -> tuple[float, float]:
        x_fraction = (
            0.5
            if max_longitude <= min_longitude
            else (longitude - min_longitude) / (max_longitude - min_longitude)
        )
        y_fraction = (
            0.5
            if max_latitude <= min_latitude
            else (latitude - min_latitude) / (max_latitude - min_latitude)
        )
        x = padding + x_fraction * (width - 2.0 * padding)
        y = padding + (1.0 - y_fraction) * (height - 2.0 * padding)
        return x, y

    estimate_by_node = {row.node: row for row in rows}
    branch_lines: list[str] = []
    for row in report.branch_rows:
        parent = estimate_by_node[row.parent_node]
        child = estimate_by_node[row.child_node]
        left_x, left_y = project(parent.longitude, parent.latitude)
        right_x, right_y = project(child.longitude, child.latitude)
        stroke = (
            "#b91c1c"
            if row.impossible_jump
            else "#c2410c"
            if row.outlier_jump
            else "#334155"
        )
        width_px = 3 if row.impossible_jump or row.outlier_jump else 1.75
        branch_lines.append(
            f'<line x1="{left_x:.1f}" y1="{left_y:.1f}" x2="{right_x:.1f}" y2="{right_y:.1f}" stroke="{stroke}" stroke-width="{width_px}" stroke-linecap="round" />'
        )

    point_marks: list[str] = []
    labels: list[str] = []
    for row in rows:
        x, y = project(row.longitude, row.latitude)
        if row.is_tip:
            fill = "#0f172a"
            radius = 4.5
        elif row.is_root:
            fill = "#a16207"
            radius = 5.0
        else:
            fill = "#0f766e"
            radius = 3.5
        point_marks.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius:.1f}" fill="{fill}" opacity="0.95" />'
        )
        if row.node_name is not None:
            labels.append(
                f'<text x="{x + 6:.1f}" y="{y - 6:.1f}" font-size="12" fill="#0f172a">{escape(row.node_name)}</text>'
            )

    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}">',
            '<rect width="100%" height="100%" fill="#f8fafc" />',
            '<text x="24" y="30" font-size="20" fill="#0f172a">Coordinate-Space Movement Review</text>',
            f'<text x="24" y="{height - 18}" font-size="12" fill="#475569">Longitude increases left to right. Latitude increases bottom to top. This is a coordinate-space view, not a projected map.</text>',
            *branch_lines,
            *point_marks,
            *labels,
            "</svg>",
        ]
    )


def _build_coordinate_html(
    report: PhylogeographicCoordinateReport,
    svg_markup: str,
) -> str:
    summary = report.summary
    return "\n".join(
        [
            "<!doctype html>",
            '<html lang="en">',
            "<head>",
            '<meta charset="utf-8">',
            "<title>Continuous phylogeography review</title>",
            "<style>",
            "body { font-family: ui-sans-serif, system-ui, sans-serif; margin: 24px; color: #0f172a; background: #f8fafc; }",
            ".metrics { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; margin-bottom: 18px; }",
            ".card { background: white; border: 1px solid #cbd5e1; border-radius: 10px; padding: 12px; }",
            ".label { font-size: 12px; color: #475569; }",
            ".value { font-size: 18px; font-weight: 700; }",
            "svg { background: white; border: 1px solid #cbd5e1; border-radius: 10px; }",
            "</style>",
            "</head>",
            "<body>",
            "<h1>Continuous Phylogeography Review</h1>",
            '<div class="metrics">',
            f'<div class="card"><div class="label">Model</div><div class="value">{escape(summary.model)}</div></div>',
            f'<div class="card"><div class="label">Analyzed Taxa</div><div class="value">{summary.analyzed_taxon_count}</div></div>',
            f'<div class="card"><div class="label">Flagged Branches</div><div class="value">{summary.flagged_branch_count}</div></div>',
            f'<div class="card"><div class="label">Maximum Jump (km)</div><div class="value">{summary.maximum_jump_km}</div></div>',
            "</div>",
            svg_markup,
            "</body>",
            "</html>",
        ]
    )


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(_stable_float(value))
