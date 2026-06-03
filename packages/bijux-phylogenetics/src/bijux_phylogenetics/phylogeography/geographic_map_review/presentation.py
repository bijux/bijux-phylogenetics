from __future__ import annotations

from html import escape
from pathlib import Path

from bijux_phylogenetics.ancestral.common import stable_value
from bijux_phylogenetics.phylogeography.region_styles import (
    build_geographic_state_color_map,
    geographic_transition_support_colors,
)

from .contracts import GeographicMapArtifact, GeographicMapReport
from .shared import format_optional_float


def render_geographic_map_html(
    report: GeographicMapReport,
    *,
    out_path: Path,
    width: int = 1200,
    height: int = 720,
    state_colors: dict[str, str] | None = None,
) -> GeographicMapArtifact:
    """Render one self-contained HTML geographic map review artifact."""
    if out_path.suffix.lower() != ".html":
        raise ValueError("geographic map output must end in .html")
    html = _build_map_html(
        report,
        width=width,
        height=height,
        state_colors=state_colors,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    return GeographicMapArtifact(
        output_path=out_path,
        format="html",
        width=width,
        height=height,
        visible_line_count=sum(row.visible for row in report.line_rows),
    )


def _build_map_html(
    report: GeographicMapReport,
    *,
    width: int,
    height: int,
    state_colors: dict[str, str] | None = None,
) -> str:
    summary = report.summary
    state_colors = state_colors or build_geographic_state_color_map(
        row.state_label for row in report.marker_rows
    )
    transition_colors = geographic_transition_support_colors()
    svg = _build_map_svg(
        report,
        width=width,
        height=height,
        state_colors=state_colors,
        transition_colors=transition_colors,
    )
    warnings_markup = "".join(
        f"<li>{escape(warning)}</li>" for warning in report.warnings
    )
    title = (
        "Continuous Geographic Map Review"
        if report.mode == "continuous"
        else "Regional Transition Map Review"
    )
    filter_note = (
        f"{format_optional_float(summary.minimum_midpoint_depth_filter)} to "
        f"{format_optional_float(summary.maximum_midpoint_depth_filter)}"
        if summary.time_filter_applied
        else "none"
    )
    return "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            f"  <title>{escape(title)}</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f7f3 0%, #e4ece8 100%); color: #163127; }",
            "    main { max-width: 1320px; margin: 0 auto; padding: 24px; }",
            "    h1 { margin: 0 0 8px; font-size: 32px; }",
            "    p { margin: 0 0 12px; line-height: 1.5; }",
            "    .grid { display: grid; grid-template-columns: 340px 1fr; gap: 20px; align-items: start; }",
            "    .panel { background: rgba(255,255,255,0.78); border: 1px solid rgba(22,49,39,0.18); border-radius: 18px; padding: 18px; box-shadow: 0 18px 44px rgba(22,49,39,0.08); }",
            "    .metric { display: flex; justify-content: space-between; gap: 12px; padding: 6px 0; border-bottom: 1px solid rgba(22,49,39,0.08); }",
            "    .metric:last-child { border-bottom: none; }",
            "    .label { color: #496559; }",
            "    ul { margin: 8px 0 0 18px; padding: 0; }",
            "    .legend-block { margin-top: 16px; }",
            "    .legend-grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 10px; margin-top: 8px; }",
            "    .legend-row { display: flex; align-items: center; gap: 8px; font-size: 13px; }",
            "    .legend-swatch { width: 14px; height: 14px; border-radius: 999px; border: 1px solid rgba(22,49,39,0.22); flex: 0 0 auto; }",
            "    .legend-line { width: 24px; height: 0; border-top-width: 3px; border-top-style: solid; flex: 0 0 auto; }",
            "    .map-shell { overflow: auto; }",
            "    svg { width: 100%; height: auto; display: block; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            f"<h1>{escape(title)}</h1>",
            "<p>Self-contained HTML map review over owned phylogenetic geography outputs. The world view uses a fixed latitude/longitude extent so spatial interpretation remains map-based rather than coordinate-space-only.</p>",
            '<div class="grid">',
            '<section class="panel">',
            '<div class="metric"><span class="label">mode</span><strong>'
            + escape(summary.mode)
            + "</strong></div>",
            '<div class="metric"><span class="label">model</span><strong>'
            + escape(summary.model)
            + "</strong></div>",
            '<div class="metric"><span class="label">tip markers</span><strong>'
            + str(summary.tip_marker_count)
            + "</strong></div>",
            '<div class="metric"><span class="label">ancestral markers</span><strong>'
            + str(summary.internal_marker_count + summary.root_marker_count)
            + "</strong></div>",
            '<div class="metric"><span class="label">visible lines</span><strong>'
            + str(summary.visible_line_count)
            + "</strong></div>",
            '<div class="metric"><span class="label">all lines</span><strong>'
            + str(summary.line_count)
            + "</strong></div>",
            '<div class="metric"><span class="label">tree depth</span><strong>'
            + str(summary.tree_depth)
            + "</strong></div>",
            '<div class="metric"><span class="label">midpoint-depth filter</span><strong>'
            + escape(filter_note)
            + "</strong></div>",
            '<div class="metric"><span class="label">excluded records</span><strong>'
            + str(summary.excluded_record_count)
            + "</strong></div>",
            '<div class="metric"><span class="label">warnings</span><strong>'
            + str(summary.warning_count)
            + "</strong></div>",
            "<p><strong>Legend.</strong> Geographic region colors are stable across the tree and map. Marker size still distinguishes tips, ancestral nodes, and the root. Transition line color follows support class, while pale dashed lines remain outside the active depth window.</p>",
            "<p><strong>Time filtering.</strong> When branch lengths represent dated depth, midpoint-depth filters can restrict the visible transition layer without refitting the reconstruction.</p>",
            '<section class="legend-block">',
            "  <strong>Region colors</strong>",
            _legend_grid(
                [
                    _legend_row_html(
                        swatch_color=color,
                        label=label,
                    )
                    for label, color in state_colors.items()
                ]
            ),
            "</section>",
            '<section class="legend-block">',
            "  <strong>Transition support</strong>",
            _legend_grid(
                [
                    _legend_row_html(
                        line_color=color,
                        label=f"{label} support",
                    )
                    for label, color in transition_colors.items()
                ]
            ),
            "</section>",
            "<ul>" + warnings_markup + "</ul>",
            "</section>",
            '<section class="panel map-shell">',
            svg,
            "</section>",
            "</div>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )


def _build_map_svg(
    report: GeographicMapReport,
    *,
    width: int,
    height: int,
    state_colors: dict[str, str],
    transition_colors: dict[str, str],
) -> str:
    margin_left = 70.0
    margin_right = 28.0
    margin_top = 28.0
    margin_bottom = 52.0
    inner_width = width - margin_left - margin_right
    inner_height = height - margin_top - margin_bottom

    def project(longitude: float, latitude: float) -> tuple[float, float]:
        x = margin_left + ((longitude + 180.0) / 360.0) * inner_width
        y = margin_top + ((90.0 - latitude) / 180.0) * inner_height
        return stable_value(x), stable_value(y)

    grid_lines: list[str] = []
    for longitude in range(-180, 181, 30):
        x, _ = project(float(longitude), 0.0)
        grid_lines.append(
            f'<line x1="{x}" y1="{margin_top}" x2="{x}" y2="{margin_top + inner_height}" stroke="#b9cfc6" stroke-width="1" stroke-dasharray="4 6" />'
        )
        grid_lines.append(
            f'<text x="{x}" y="{height - 16}" font-size="11" text-anchor="middle" fill="#557567">{longitude}</text>'
        )
    for latitude in range(-60, 91, 30):
        _, y = project(0.0, float(latitude))
        grid_lines.append(
            f'<line x1="{margin_left}" y1="{y}" x2="{margin_left + inner_width}" y2="{y}" stroke="#b9cfc6" stroke-width="1" stroke-dasharray="4 6" />'
        )
        grid_lines.append(
            f'<text x="18" y="{y + 4}" font-size="11" fill="#557567">{latitude}</text>'
        )

    line_elements: list[str] = []
    for row in report.line_rows:
        left_x, left_y = project(row.source_longitude, row.source_latitude)
        right_x, right_y = project(row.target_longitude, row.target_latitude)
        confidence_class = row.flag_codes[0] if row.flag_codes else ""
        stroke = (
            "#0c5a74"
            if row.line_kind == "continuous_movement"
            else transition_colors.get(confidence_class, "#186b3d")
        )
        opacity = "0.88" if row.visible else "0.22"
        dasharray = "" if row.visible else "6 7"
        line_elements.append(
            "\n".join(
                [
                    f'<line x1="{left_x}" y1="{left_y}" x2="{right_x}" y2="{right_y}" stroke="{stroke}" stroke-width="{3.6 if row.visible else 2.0}" opacity="{opacity}"'
                    + (f' stroke-dasharray="{dasharray}"' if dasharray else "")
                    + " >",
                    "<title>"
                    + escape(
                        f"{row.source_label} -> {row.target_label}; distance={row.distance_km} km; midpoint_depth={format_optional_float(row.midpoint_depth)}; support={row.support}; transition={row.state_transition or 'continuous'}"
                    )
                    + "</title>",
                    "</line>",
                ]
            )
        )

    marker_elements: list[str] = []
    for row in report.marker_rows:
        x, y = project(row.longitude, row.latitude)
        fill = state_colors.get(row.state_label, "#0f6090")
        radius = 5.0 if row.is_tip else 7.0
        stroke = "#ffffff"
        stroke_width = 1.4
        if row.is_root:
            radius = 8.0
            stroke = "#17261f"
            stroke_width = 2.1
        opacity = 0.95 if row.active_line_count or row.is_tip else 0.78
        marker_elements.append(
            "\n".join(
                [
                    f'<circle cx="{x}" cy="{y}" r="{radius}" fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}" opacity="{opacity}">',
                    "<title>"
                    + escape(
                        f"{row.label}; lat={row.latitude}; lon={row.longitude}; state={row.state_label or 'continuous'}; confidence={row.confidence}; uncertainty_km={format_optional_float(row.uncertainty_km)}"
                    )
                    + "</title>",
                    "</circle>",
                ]
            )
        )
    return "\n".join(
        [
            f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Geographic reconstruction map">',
            f'<rect x="{margin_left}" y="{margin_top}" width="{inner_width}" height="{inner_height}" rx="24" fill="#edf5f2" stroke="#7aa28f" stroke-width="1.5" />',
            *grid_lines,
            *line_elements,
            *marker_elements,
            f'<text x="{margin_left}" y="18" font-size="13" fill="#355648">longitude / latitude map extent</text>',
            "</svg>",
        ]
    )


def _legend_row_html(
    *,
    label: str,
    swatch_color: str | None = None,
    line_color: str | None = None,
) -> str:
    icon = ""
    if swatch_color is not None:
        icon = f'<span class="legend-swatch" style="background:{escape(swatch_color)};"></span>'
    elif line_color is not None:
        icon = f'<span class="legend-line" style="border-top-color:{escape(line_color)};"></span>'
    return f'<div class="legend-row">{icon}<span>{escape(label)}</span></div>'


def _legend_grid(rows: list[str]) -> str:
    if not rows:
        return "<p>No legend entries.</p>"
    return '<div class="legend-grid">' + "".join(rows) + "</div>"
