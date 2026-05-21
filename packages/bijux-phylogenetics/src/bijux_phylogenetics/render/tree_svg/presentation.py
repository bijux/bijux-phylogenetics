"""SVG document assembly for rendered tree figures."""

from __future__ import annotations

from html import escape

from .contracts import AnnotationStrip
from .render_state import TreeSvgRenderState
from .shared import continuous_color, nice_scale_bar_length
from .support_audit import format_branch_value


def build_tree_svg_document(
    *,
    width: int,
    height: int,
    layout: str,
    max_distance: float,
    tree_width: float,
    left_margin: int,
    top_margin: int,
    metadata_strips: list[AnnotationStrip],
    heatmap_columns: list[AnnotationStrip],
    categorical_color_map: dict[str, str],
    continuous_values: list[float],
    continuous_min: float,
    continuous_max: float,
    state: TreeSvgRenderState,
) -> tuple[str, bool, float | None]:
    """Build the final SVG document and reviewer-facing legend surfaces."""

    scale_bar = ""
    has_scale_bar = layout == "phylogram" and max_distance > 0
    scale_bar_length: float | None = None
    if has_scale_bar:
        scale_length = nice_scale_bar_length(max_distance)
        scale_bar_length = scale_length
        scale_start = left_margin
        scale_end = left_margin + (scale_length / max_distance) * tree_width
        scale_y = height - 28
        scale_bar = (
            f'<line x1="{scale_start:.1f}" y1="{scale_y:.1f}" x2="{scale_end:.1f}" y2="{scale_y:.1f}" class="scale-bar"/>'
            f'<line x1="{scale_start:.1f}" y1="{scale_y - 6:.1f}" x2="{scale_start:.1f}" y2="{scale_y + 6:.1f}" class="scale-bar"/>'
            f'<line x1="{scale_end:.1f}" y1="{scale_y - 6:.1f}" x2="{scale_end:.1f}" y2="{scale_y + 6:.1f}" class="scale-bar"/>'
            f'<text x="{(scale_start + scale_end) / 2:.1f}" y="{scale_y - 10:.1f}" class="scale-label">{escape(format_branch_value(scale_length))}</text>'
        )

    categorical_legend = ""
    if categorical_color_map:
        legend_items = []
        legend_x = width - 220
        legend_y = 44
        for index, category in enumerate(sorted(categorical_color_map)):
            item_y = legend_y + 24 + index * 22
            legend_items.append(
                f'<circle cx="{legend_x:.1f}" cy="{item_y - 4:.1f}" r="6" fill="{categorical_color_map[category]}" class="trait-marker"/>'
                f'<text x="{legend_x + 16:.1f}" y="{item_y:.1f}" class="legend-label">{escape(category)}</text>'
            )
        categorical_legend = (
            f'<text x="{legend_x - 2:.1f}" y="{legend_y:.1f}" class="legend-title">categorical trait</text>'
            + "".join(legend_items)
        )

    continuous_legend = ""
    if continuous_values:
        legend_x = width - 220
        legend_y = height - 60
        continuous_legend = (
            f'<text x="{legend_x - 2:.1f}" y="{legend_y - 10:.1f}" class="legend-title">continuous trait</text>'
            f'<defs><linearGradient id="continuous-trait-gradient" x1="0%" y1="0%" x2="100%" y2="0%">'
            f'<stop offset="0%" stop-color="{continuous_color(continuous_min, continuous_min, continuous_max)}"/>'
            f'<stop offset="100%" stop-color="{continuous_color(continuous_max, continuous_min, continuous_max)}"/>'
            f"</linearGradient></defs>"
            f'<rect x="{legend_x:.1f}" y="{legend_y:.1f}" width="120" height="12" rx="6" fill="url(#continuous-trait-gradient)"/>'
            f'<text x="{legend_x:.1f}" y="{legend_y + 28:.1f}" class="legend-label">{escape(format_branch_value(continuous_min))}</text>'
            f'<text x="{legend_x + 88:.1f}" y="{legend_y + 28:.1f}" class="legend-label">{escape(format_branch_value(continuous_max))}</text>'
        )

    metadata_strip_headers = ""
    if metadata_strips and layout != "circular":
        header_fragments = []
        base_x = left_margin + tree_width + 298
        header_y = top_margin - 10
        for strip_index, strip in enumerate(metadata_strips):
            header_x = base_x + strip_index * 24
            header_fragments.append(
                f'<text x="{header_x:.1f}" y="{header_y:.1f}" transform="rotate(-35 {header_x:.1f} {header_y:.1f})" class="strip-header">{escape(strip.name)}</text>'
            )
        metadata_strip_headers = "".join(header_fragments)

    heatmap_headers = ""
    if heatmap_columns and layout != "circular":
        header_fragments = []
        base_x = left_margin + tree_width + 326 + len(metadata_strips) * 24
        header_y = top_margin - 10
        for column_index, column in enumerate(heatmap_columns):
            header_x = base_x + column_index * 22
            header_fragments.append(
                f'<text x="{header_x:.1f}" y="{header_y:.1f}" transform="rotate(-35 {header_x:.1f} {header_y:.1f})" class="strip-header">{escape(column.name)}</text>'
            )
        heatmap_headers = "".join(header_fragments)

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="phylogenetic tree">
  <style>
    .panel {{ fill: #f7fbfa; stroke: #d7e3e1; stroke-width: 1; rx: 18; ry: 18; }}
    .branch {{ stroke: #0f172a; stroke-width: 2.2; stroke-linecap: round; fill: none; }}
    .branch-colored {{ stroke-width: 2.8; }}
    .scale-bar {{ stroke: #0f172a; stroke-width: 2; stroke-linecap: round; }}
    .collapsed-clade {{ fill-opacity: 0.95; }}
    .metadata-strip-cell {{ stroke: #f8fafc; stroke-width: 1; }}
    .heatmap-cell {{ stroke: #f8fafc; stroke-width: 1; }}
    .internal-pie-slice {{ stroke: #f8fafc; stroke-width: 0.9; }}
    .tip-label {{ fill: #0f172a; font: 16px "Avenir Next", "Segoe UI", sans-serif; }}
    .support-label {{ fill: #0f766e; font: 12px "Avenir Next", "Segoe UI", sans-serif; }}
    .internal-annotation-label {{ fill: #6d28d9; font: 12px "Avenir Next", "Segoe UI", sans-serif; }}
    .legend-title {{ fill: #334155; font: 12px "Avenir Next", "Segoe UI", sans-serif; text-transform: uppercase; letter-spacing: 0.08em; }}
    .legend-label {{ fill: #334155; font: 13px "Avenir Next", "Segoe UI", sans-serif; }}
    .scale-label {{ fill: #334155; text-anchor: middle; font: 13px "Avenir Next", "Segoe UI", sans-serif; }}
    .strip-header {{ fill: #475569; font: 11px "Avenir Next", "Segoe UI", sans-serif; }}
  </style>
  <rect x="1" y="1" width="{width - 2}" height="{height - 2}" class="panel" />
  {"".join(state.lines)}
  {"".join(state.texts)}
  {"".join(state.overlays)}
  {scale_bar}
  {categorical_legend}
  {continuous_legend}
  {metadata_strip_headers}
  {heatmap_headers}
</svg>
"""
    return svg, has_scale_bar, scale_bar_length
