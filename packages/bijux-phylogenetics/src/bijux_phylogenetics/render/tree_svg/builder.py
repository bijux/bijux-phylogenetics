"""Build deterministic SVG tree figures from one phylogenetic tree input."""

from __future__ import annotations

from html import escape
import math
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import _load_tree

from .contracts import AnnotationStrip, TreeRenderResult
from .rectangular_layout import render_rectangular_tree_layout
from .render_state import TreeSvgRenderState
from .shared import (
    categorical_color_map as build_categorical_color_map,
    continuous_color,
    count_subtree_leaves,
    count_visible_leaves,
    is_collapsed_node,
    is_numeric_strings,
    max_visible_depth,
    max_visible_distance,
    nice_scale_bar_length,
    node_signature,
    polar_point,
    svg_pie_slices,
)
from .support_audit import coerce_support_label, format_branch_value


def render_tree_svg(
    tree_path: Path,
    *,
    out_path: Path,
    labels: dict[str, str] | None = None,
    layout: str = "cladogram",
    show_support_values: bool = False,
    categorical_traits: dict[str, str] | None = None,
    continuous_traits: dict[str, float] | None = None,
    metadata_strips: list[AnnotationStrip] | None = None,
    heatmap_columns: list[AnnotationStrip] | None = None,
    collapsed_clades: list[str] | None = None,
    internal_annotations: dict[str, str] | None = None,
    internal_annotation_colors: dict[str, str] | None = None,
    branch_colors: dict[str, str] | None = None,
    internal_pies: dict[str, dict[str, float]] | None = None,
    internal_pie_colors: dict[str, str] | None = None,
    validated_support_labels: dict[str, str] | None = None,
    support_validation_warnings: list[str] | None = None,
) -> TreeRenderResult:
    """Render a deterministic SVG tree as a cladogram, phylogram, or circular tree."""
    if layout not in {"cladogram", "phylogram", "circular"}:
        raise ValueError(f"unsupported tree layout: {layout}")

    tree = _load_tree(tree_path)
    labels = labels or {}
    collapsed_clade_names = {name for name in (collapsed_clades or []) if name}
    row_height = 56
    left_margin = 48
    right_margin = 560
    top_margin = 40
    bottom_margin = 72
    horizontal_step = 150
    scale_width = 520
    leaf_count = count_visible_leaves(tree.root, collapsed_clade_names)
    max_depth = max(max_visible_depth(tree.root, collapsed_clade_names), 1)
    max_distance = (
        max_visible_distance(tree.root, collapsed_clade_names, 0.0)
        if layout in {"phylogram", "circular"}
        else 0.0
    )

    if layout == "circular":
        tree_radius = 320
        width = tree_radius * 2 + 280
        height = tree_radius * 2 + 120
        tree_width = tree_radius
    elif layout == "phylogram" and max_distance > 0:
        tree_width = scale_width
        width = left_margin + tree_width + right_margin
        height = top_margin + bottom_margin + row_height * max(leaf_count, 1)
    else:
        tree_width = horizontal_step * (max_depth + 1)
        width = left_margin + tree_width + right_margin
        height = top_margin + bottom_margin + row_height * max(leaf_count, 1)

    state = TreeSvgRenderState()
    validated_support_labels = validated_support_labels or {}
    support_validation_warnings = support_validation_warnings or []
    categorical_traits = categorical_traits or {}
    categorical_color_map = build_categorical_color_map(categorical_traits)
    continuous_traits = continuous_traits or {}
    continuous_values = list(continuous_traits.values())
    continuous_min = min(continuous_values) if continuous_values else 0.0
    continuous_max = max(continuous_values) if continuous_values else 0.0
    metadata_strips = metadata_strips or []
    metadata_strip_colors = [
        build_categorical_color_map(strip.values) for strip in metadata_strips
    ]
    heatmap_columns = heatmap_columns or []
    internal_annotations = internal_annotations or {}
    internal_annotation_colors = internal_annotation_colors or {}
    branch_colors = branch_colors or {}
    internal_pies = internal_pies or {}
    internal_pie_colors = internal_pie_colors or {}
    heatmap_specs: list[tuple[str, dict[str, str], float, float]] = []
    for column in heatmap_columns:
        observed_values = [value for value in column.values.values() if value]
        if is_numeric_strings(observed_values):
            numeric_values = [float(value) for value in observed_values]
            heatmap_specs.append(
                ("numeric", {}, min(numeric_values), max(numeric_values))
            )
        else:
            heatmap_specs.append(
                ("categorical", build_categorical_color_map(column.values), 0.0, 0.0)
            )

    def visit_circular() -> None:
        center_x = width / 2
        center_y = height / 2
        radius = min(width, height) / 2 - 80
        angle_cache: dict[str, float] = {}

        def radial_distance(depth: int, distance: float) -> float:
            if max_distance > 0:
                return (distance / max_distance) * radius
            return (depth / max(max_depth, 1)) * radius

        def assign_angles(node: TreeNode) -> tuple[float, float]:
            if node.is_leaf() or is_collapsed_node(node, collapsed_clade_names):
                angle = (
                    2 * math.pi * state.next_leaf_index / max(leaf_count, 1)
                ) - math.pi / 2
                state.next_leaf_index += 1
                angle_cache[node.node_id or node_signature(node)] = angle
                return angle, angle
            ranges = [assign_angles(child) for child in node.children]
            start_angle = min(start for start, _ in ranges)
            end_angle = max(end for _, end in ranges)
            angle_cache[node.node_id or node_signature(node)] = (
                start_angle + end_angle
            ) / 2
            return start_angle, end_angle

        def draw(node: TreeNode, depth: int, distance: float) -> tuple[float, float]:
            branch_distance = distance + float(node.branch_length or 0.0)
            radial = radial_distance(
                depth, branch_distance if node is not tree.root else distance
            )
            if node.is_leaf() or is_collapsed_node(node, collapsed_clade_names):
                angle = angle_cache[node.node_id or node_signature(node)]
                label = labels.get(node.name or "", node.name or "")
                if is_collapsed_node(node, collapsed_clade_names):
                    label = (
                        f"{node.name or 'collapsed clade'}"
                        f" ({count_subtree_leaves(node)} tips)"
                    )
                    state.rendered_collapsed_clades += 1
                elif node.name and node.name not in labels and labels:
                    state.missing_labels.append(node.name)
                anchor = "start" if math.cos(angle) >= 0 else "end"
                label_point = polar_point(center_x, center_y, radial + 18, angle)
                state.texts.append(
                    f'<text x="{label_point.x:.1f}" y="{label_point.y + 5:.1f}" text-anchor="{anchor}" class="tip-label">{escape(label)}</text>'
                )
                if is_collapsed_node(node, collapsed_clade_names):
                    node_point = polar_point(center_x, center_y, radial, angle)
                    left_point = polar_point(
                        center_x, center_y, radial + 18, angle - 0.08
                    )
                    right_point = polar_point(
                        center_x, center_y, radial + 18, angle + 0.08
                    )
                    state.overlays.append(
                        f'<polygon points="{node_point.x:.1f},{node_point.y:.1f} {left_point.x:.1f},{left_point.y:.1f} {right_point.x:.1f},{right_point.y:.1f}" fill="#cbd5e1" stroke="#475569" stroke-width="1.5" class="collapsed-clade"/>'
                    )
                trait_value = categorical_traits.get(node.name or "")
                if trait_value:
                    marker_point = polar_point(center_x, center_y, radial + 10, angle)
                    state.overlays.append(
                        f'<circle cx="{marker_point.x:.1f}" cy="{marker_point.y:.1f}" r="6" fill="{categorical_color_map[trait_value]}" class="trait-marker"/>'
                    )
                    state.rendered_categorical_trait_count += 1
                continuous_value = continuous_traits.get(node.name or "")
                if continuous_value is not None:
                    marker_point = polar_point(center_x, center_y, radial + 10, angle)
                    state.overlays.append(
                        f'<circle cx="{marker_point.x:.1f}" cy="{marker_point.y:.1f}" r="6" fill="{continuous_color(continuous_value, continuous_min, continuous_max)}" class="trait-marker"/>'
                    )
                    state.rendered_continuous_trait_count += 1
                return angle, radial

            child_positions = [
                draw(child, depth + 1, branch_distance) for child in node.children
            ]
            start_angle = min(angle for angle, _ in child_positions)
            end_angle = max(angle for angle, _ in child_positions)
            if radial > 0 and start_angle != end_angle:
                arc_start = polar_point(center_x, center_y, radial, start_angle)
                arc_end = polar_point(center_x, center_y, radial, end_angle)
                large_arc = 1 if end_angle - start_angle > math.pi else 0
                state.lines.append(
                    f'<path d="M {arc_start.x:.1f} {arc_start.y:.1f} A {radial:.1f} {radial:.1f} 0 {large_arc} 1 {arc_end.x:.1f} {arc_end.y:.1f}" class="branch"/>'
                )
            for child in node.children:
                child_angle = angle_cache[child.node_id or node_signature(child)]
                child_branch_distance = branch_distance + float(
                    child.branch_length or 0.0
                )
                child_radial = radial_distance(depth + 1, child_branch_distance)
                radial_start = polar_point(center_x, center_y, radial, child_angle)
                radial_end = polar_point(center_x, center_y, child_radial, child_angle)
                branch_signature = node_signature(child)
                branch_color = branch_colors.get(branch_signature)
                branch_attributes = (
                    f' stroke="{branch_color}" class="branch branch-colored"'
                    if branch_color is not None
                    else ' class="branch"'
                )
                state.lines.append(
                    f'<line x1="{radial_start.x:.1f}" y1="{radial_start.y:.1f}" x2="{radial_end.x:.1f}" y2="{radial_end.y:.1f}"{branch_attributes}/>'
                )
                if branch_color is not None:
                    state.rendered_branch_color_count += 1
            support_label = None
            if show_support_values:
                support_label = validated_support_labels.get(node_signature(node))
                if support_label is None:
                    support_label = coerce_support_label(node.name)
            if support_label is not None and node is not tree.root:
                node_angle = angle_cache[node.node_id or node_signature(node)]
                support_point = polar_point(
                    center_x, center_y, radial + 14, node_angle
                )
                state.texts.append(
                    f'<text x="{support_point.x:.1f}" y="{support_point.y + 4:.1f}" class="support-label">{escape(support_label)}</text>'
                )
                state.rendered_support_count += 1
            annotation = internal_annotations.get(node_signature(node))
            pie_segments = internal_pies.get(node_signature(node))
            if pie_segments and not node.is_leaf():
                node_point = polar_point(
                    center_x,
                    center_y,
                    radial,
                    angle_cache[node.node_id or node_signature(node)],
                )
                state.overlays.extend(
                    svg_pie_slices(
                        node_point.x,
                        node_point.y,
                        7.0,
                        pie_segments,
                        internal_pie_colors,
                    )
                )
                state.rendered_internal_pie_count += 1
            if annotation and not node.is_leaf():
                node_angle = angle_cache[node.node_id or node_signature(node)]
                annotation_point = polar_point(
                    center_x, center_y, radial + 12, node_angle
                )
                annotation_color = internal_annotation_colors.get(
                    node_signature(node), "#7c3aed"
                )
                state.overlays.append(
                    f'<circle cx="{annotation_point.x:.1f}" cy="{annotation_point.y:.1f}" r="4.5" fill="{annotation_color}" class="internal-annotation-dot"/>'
                )
                state.texts.append(
                    f'<text x="{annotation_point.x + 10:.1f}" y="{annotation_point.y + 4:.1f}" class="internal-annotation-label">{escape(annotation)}</text>'
                )
                state.rendered_internal_annotation_count += 1
            return angle_cache[node.node_id or node_signature(node)], radial

        state.next_leaf_index = 0
        assign_angles(tree.root)
        draw(tree.root, 0, 0.0)

    if layout == "circular":
        visit_circular()
    else:
        render_rectangular_tree_layout(
            tree,
            layout=layout,
            state=state,
            labels=labels,
            collapsed_clade_names=collapsed_clade_names,
            row_height=row_height,
            top_margin=top_margin,
            left_margin=left_margin,
            horizontal_step=horizontal_step,
            tree_width=tree_width,
            max_distance=max_distance,
            categorical_traits=categorical_traits,
            categorical_color_map=categorical_color_map,
            continuous_traits=continuous_traits,
            continuous_min=continuous_min,
            continuous_max=continuous_max,
            metadata_strips=metadata_strips,
            metadata_strip_colors=metadata_strip_colors,
            heatmap_columns=heatmap_columns,
            heatmap_specs=heatmap_specs,
            internal_annotations=internal_annotations,
            internal_annotation_colors=internal_annotation_colors,
            branch_colors=branch_colors,
            internal_pies=internal_pies,
            internal_pie_colors=internal_pie_colors,
            show_support_values=show_support_values,
            validated_support_labels=validated_support_labels,
        )

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
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(svg, encoding="utf-8")
    return TreeRenderResult(
        output_path=out_path,
        format="svg",
        layout=layout,
        tip_count=tree.tip_count,
        visible_tip_count=leaf_count,
        label_count=len(state.texts),
        has_scale_bar=has_scale_bar,
        scale_bar_length=scale_bar_length,
        max_branch_distance=max_distance
        if layout in {"phylogram", "circular"}
        else None,
        rendered_support_count=state.rendered_support_count,
        support_labels_validated=bool(validated_support_labels),
        support_validation_warnings=list(support_validation_warnings),
        rendered_categorical_trait_count=state.rendered_categorical_trait_count,
        rendered_continuous_trait_count=state.rendered_continuous_trait_count,
        rendered_metadata_strip_count=len(metadata_strips),
        rendered_heatmap_column_count=len(heatmap_columns),
        rendered_internal_annotation_count=state.rendered_internal_annotation_count,
        rendered_branch_color_count=state.rendered_branch_color_count,
        rendered_internal_pie_count=state.rendered_internal_pie_count,
        collapsed_clade_count=state.rendered_collapsed_clades,
        missing_metadata_labels=sorted(set(state.missing_labels)),
    )
