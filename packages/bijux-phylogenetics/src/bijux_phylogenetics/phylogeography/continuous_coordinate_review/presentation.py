from __future__ import annotations

from html import escape
from pathlib import Path

from bijux_phylogenetics.phylogeography.continuous_coordinate_review.contracts import (
    CoordinateMovementVisualization,
    PhylogeographicCoordinateReport,
)


def render_coordinate_movement_visualization(
    report: PhylogeographicCoordinateReport,
    *,
    out_path: Path,
    width: int = 960,
    height: int = 640,
) -> CoordinateMovementVisualization:
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
