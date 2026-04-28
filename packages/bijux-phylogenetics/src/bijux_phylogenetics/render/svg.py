from __future__ import annotations

from dataclasses import dataclass
from html import escape
import math
from pathlib import Path

from bijux_phylogenetics.core.tree import TreeNode
from bijux_phylogenetics.diagnostics.validation import _load_tree


@dataclass(slots=True)
class TreeRenderResult:
    output_path: Path
    format: str
    layout: str
    tip_count: int
    label_count: int
    has_scale_bar: bool
    rendered_support_count: int
    rendered_categorical_trait_count: int
    rendered_continuous_trait_count: int
    rendered_metadata_strip_count: int
    rendered_heatmap_column_count: int
    collapsed_clade_count: int
    missing_metadata_labels: list[str]


@dataclass(frozen=True, slots=True)
class _Point:
    x: float
    y: float


_CATEGORICAL_PALETTE = (
    "#0f766e",
    "#1d4ed8",
    "#c2410c",
    "#7c3aed",
    "#b91c1c",
    "#047857",
    "#a16207",
    "#0f172a",
)


def _count_render_leaves(node: TreeNode) -> int:
    if node.is_leaf():
        return 1
    return sum(_count_render_leaves(child) for child in node.children)


def _max_depth(node: TreeNode, depth: int = 0) -> int:
    if node.is_leaf():
        return depth
    return max(_max_depth(child, depth + 1) for child in node.children)


def _max_distance(node: TreeNode, distance: float = 0.0) -> float:
    next_distance = distance + float(node.branch_length or 0.0)
    if node.is_leaf():
        return next_distance
    return max(_max_distance(child, next_distance) for child in node.children)


def _nice_scale_bar_length(max_distance: float) -> float:
    if max_distance <= 0:
        return 0.0
    exponent = math.floor(math.log10(max_distance))
    base = 10**exponent
    for factor in (1.0, 0.5, 0.2, 0.1):
        candidate = base * factor
        if candidate <= max_distance / 3:
            return candidate
    return base / 10


def _format_branch_value(value: float) -> str:
    return format(round(value, 15), ".15g")


def _coerce_support_label(raw: str | None) -> str | None:
    if raw is None or not raw.strip():
        return None
    try:
        return _format_branch_value(float(raw))
    except ValueError:
        return None


def _polar_point(center_x: float, center_y: float, radius: float, angle_radians: float) -> _Point:
    return _Point(
        x=center_x + radius * math.cos(angle_radians),
        y=center_y + radius * math.sin(angle_radians),
    )


def _categorical_color_map(values: dict[str, str]) -> dict[str, str]:
    categories = sorted({value for value in values.values() if value})
    return {
        category: _CATEGORICAL_PALETTE[index % len(_CATEGORICAL_PALETTE)]
        for index, category in enumerate(categories)
    }


def render_tree_svg(
    tree_path: Path,
    *,
    out_path: Path,
    labels: dict[str, str] | None = None,
    layout: str = "cladogram",
    show_support_values: bool = False,
    categorical_traits: dict[str, str] | None = None,
) -> TreeRenderResult:
    """Render a deterministic SVG tree as a cladogram, phylogram, or circular tree."""
    if layout not in {"cladogram", "phylogram", "circular"}:
        raise ValueError(f"unsupported tree layout: {layout}")

    tree = _load_tree(tree_path)
    labels = labels or {}
    row_height = 56
    left_margin = 48
    right_margin = 320
    top_margin = 40
    bottom_margin = 72
    horizontal_step = 150
    scale_width = 520
    leaf_count = _count_render_leaves(tree.root)
    max_depth = max(_max_depth(tree.root), 1)
    max_distance = _max_distance(tree.root, 0.0) if layout in {"phylogram", "circular"} else 0.0

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

    lines: list[str] = []
    texts: list[str] = []
    overlays: list[str] = []
    missing_labels: list[str] = []
    next_leaf_index = 0
    rendered_support_count = 0
    rendered_categorical_trait_count = 0
    categorical_traits = categorical_traits or {}
    categorical_color_map = _categorical_color_map(categorical_traits)

    def node_x(depth: int, distance: float) -> float:
        if layout == "phylogram" and max_distance > 0:
            return left_margin + (distance / max_distance) * tree_width
        return left_margin + depth * horizontal_step

    def visit_rectangular(node: TreeNode, depth: int, distance: float) -> _Point:
        nonlocal next_leaf_index
        nonlocal rendered_support_count
        nonlocal rendered_categorical_trait_count
        branch_distance = distance + float(node.branch_length or 0.0)
        x = node_x(depth, branch_distance if node is not tree.root else distance)
        if node.is_leaf():
            y = top_margin + next_leaf_index * row_height + row_height / 2
            next_leaf_index += 1
            label = labels.get(node.name or "", node.name or "")
            if node.name and node.name not in labels and labels:
                missing_labels.append(node.name)
            label_x = x + 18
            texts.append(
                f'<text x="{label_x:.1f}" y="{y + 5:.1f}" class="tip-label">{escape(label)}</text>'
            )
            trait_value = categorical_traits.get(node.name or "")
            if trait_value:
                marker_x = left_margin + tree_width + 190
                overlays.append(
                    f'<circle cx="{marker_x:.1f}" cy="{y - 4:.1f}" r="7" fill="{categorical_color_map[trait_value]}" class="trait-marker"/>'
                )
                rendered_categorical_trait_count += 1
            return _Point(x=x, y=y)

        child_points = [visit_rectangular(child, depth + 1, branch_distance) for child in node.children]
        y = sum(point.y for point in child_points) / len(child_points)
        min_y = min(point.y for point in child_points)
        max_y = max(point.y for point in child_points)
        lines.append(
            f'<line x1="{x:.1f}" y1="{min_y:.1f}" x2="{x:.1f}" y2="{max_y:.1f}" class="branch"/>'
        )
        for child_point in child_points:
            lines.append(
                f'<line x1="{x:.1f}" y1="{child_point.y:.1f}" x2="{child_point.x:.1f}" y2="{child_point.y:.1f}" class="branch"/>'
            )
        support_label = _coerce_support_label(node.name) if show_support_values else None
        if support_label is not None and node is not tree.root:
            texts.append(
                f'<text x="{x + 8:.1f}" y="{y - 8:.1f}" class="support-label">{escape(support_label)}</text>'
            )
            rendered_support_count += 1
        return _Point(x=x, y=y)

    def visit_circular() -> None:
        nonlocal next_leaf_index
        nonlocal rendered_support_count
        center_x = width / 2
        center_y = height / 2
        radius = min(width, height) / 2 - 80
        angle_cache: dict[int, float] = {}

        def radial_distance(depth: int, distance: float) -> float:
            if max_distance > 0:
                return (distance / max_distance) * radius
            return (depth / max(max_depth, 1)) * radius

        def assign_angles(node: TreeNode) -> tuple[float, float]:
            nonlocal next_leaf_index
            if node.is_leaf():
                angle = (2 * math.pi * next_leaf_index / max(leaf_count, 1)) - math.pi / 2
                next_leaf_index += 1
                angle_cache[id(node)] = angle
                return angle, angle
            ranges = [assign_angles(child) for child in node.children]
            start_angle = min(start for start, _ in ranges)
            end_angle = max(end for _, end in ranges)
            angle_cache[id(node)] = (start_angle + end_angle) / 2
            return start_angle, end_angle

        def draw(node: TreeNode, depth: int, distance: float) -> tuple[float, float]:
            nonlocal rendered_support_count
            nonlocal rendered_categorical_trait_count
            branch_distance = distance + float(node.branch_length or 0.0)
            radial = radial_distance(depth, branch_distance if node is not tree.root else distance)
            if node.is_leaf():
                angle = angle_cache[id(node)]
                label = labels.get(node.name or "", node.name or "")
                if node.name and node.name not in labels and labels:
                    missing_labels.append(node.name)
                anchor = "start" if math.cos(angle) >= 0 else "end"
                label_point = _polar_point(center_x, center_y, radial + 18, angle)
                texts.append(
                    f'<text x="{label_point.x:.1f}" y="{label_point.y + 5:.1f}" text-anchor="{anchor}" class="tip-label">{escape(label)}</text>'
                )
                trait_value = categorical_traits.get(node.name or "")
                if trait_value:
                    marker_point = _polar_point(center_x, center_y, radial + 10, angle)
                    overlays.append(
                        f'<circle cx="{marker_point.x:.1f}" cy="{marker_point.y:.1f}" r="6" fill="{categorical_color_map[trait_value]}" class="trait-marker"/>'
                    )
                    rendered_categorical_trait_count += 1
                return angle, radial

            child_positions = [draw(child, depth + 1, branch_distance) for child in node.children]
            start_angle = min(angle for angle, _ in child_positions)
            end_angle = max(angle for angle, _ in child_positions)
            if radial > 0 and start_angle != end_angle:
                arc_start = _polar_point(center_x, center_y, radial, start_angle)
                arc_end = _polar_point(center_x, center_y, radial, end_angle)
                large_arc = 1 if end_angle - start_angle > math.pi else 0
                lines.append(
                    f'<path d="M {arc_start.x:.1f} {arc_start.y:.1f} A {radial:.1f} {radial:.1f} 0 {large_arc} 1 {arc_end.x:.1f} {arc_end.y:.1f}" class="branch"/>'
                )
            for child in node.children:
                child_angle = angle_cache[id(child)]
                child_branch_distance = branch_distance + float(child.branch_length or 0.0)
                child_radial = radial_distance(depth + 1, child_branch_distance)
                radial_start = _polar_point(center_x, center_y, radial, child_angle)
                radial_end = _polar_point(center_x, center_y, child_radial, child_angle)
                lines.append(
                    f'<line x1="{radial_start.x:.1f}" y1="{radial_start.y:.1f}" x2="{radial_end.x:.1f}" y2="{radial_end.y:.1f}" class="branch"/>'
                )
            support_label = _coerce_support_label(node.name) if show_support_values else None
            if support_label is not None and node is not tree.root:
                node_angle = angle_cache[id(node)]
                support_point = _polar_point(center_x, center_y, radial + 14, node_angle)
                texts.append(
                    f'<text x="{support_point.x:.1f}" y="{support_point.y + 4:.1f}" class="support-label">{escape(support_label)}</text>'
                )
                rendered_support_count += 1
            return angle_cache[id(node)], radial

        next_leaf_index = 0
        assign_angles(tree.root)
        draw(tree.root, 0, 0.0)

    if layout == "circular":
        visit_circular()
    else:
        visit_rectangular(tree.root, 0, 0.0)

    scale_bar = ""
    has_scale_bar = layout == "phylogram" and max_distance > 0
    if has_scale_bar:
        scale_length = _nice_scale_bar_length(max_distance)
        scale_start = left_margin
        scale_end = left_margin + (scale_length / max_distance) * tree_width
        scale_y = height - 28
        scale_bar = (
            f'<line x1="{scale_start:.1f}" y1="{scale_y:.1f}" x2="{scale_end:.1f}" y2="{scale_y:.1f}" class="scale-bar"/>'
            f'<line x1="{scale_start:.1f}" y1="{scale_y - 6:.1f}" x2="{scale_start:.1f}" y2="{scale_y + 6:.1f}" class="scale-bar"/>'
            f'<line x1="{scale_end:.1f}" y1="{scale_y - 6:.1f}" x2="{scale_end:.1f}" y2="{scale_y + 6:.1f}" class="scale-bar"/>'
            f'<text x="{(scale_start + scale_end) / 2:.1f}" y="{scale_y - 10:.1f}" class="scale-label">{escape(_format_branch_value(scale_length))}</text>'
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

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="phylogenetic tree">
  <style>
    .panel {{ fill: #f7fbfa; stroke: #d7e3e1; stroke-width: 1; rx: 18; ry: 18; }}
    .branch {{ stroke: #0f172a; stroke-width: 2.2; stroke-linecap: round; fill: none; }}
    .scale-bar {{ stroke: #0f172a; stroke-width: 2; stroke-linecap: round; }}
    .tip-label {{ fill: #0f172a; font: 16px "Avenir Next", "Segoe UI", sans-serif; }}
    .support-label {{ fill: #0f766e; font: 12px "Avenir Next", "Segoe UI", sans-serif; }}
    .legend-title {{ fill: #334155; font: 12px "Avenir Next", "Segoe UI", sans-serif; text-transform: uppercase; letter-spacing: 0.08em; }}
    .legend-label {{ fill: #334155; font: 13px "Avenir Next", "Segoe UI", sans-serif; }}
    .scale-label {{ fill: #334155; text-anchor: middle; font: 13px "Avenir Next", "Segoe UI", sans-serif; }}
  </style>
  <rect x="1" y="1" width="{width - 2}" height="{height - 2}" class="panel" />
  {''.join(lines)}
  {''.join(texts)}
  {''.join(overlays)}
  {scale_bar}
  {categorical_legend}
</svg>
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(svg, encoding="utf-8")
    return TreeRenderResult(
        output_path=out_path,
        format="svg",
        layout=layout,
        tip_count=tree.tip_count,
        label_count=len(texts),
        has_scale_bar=has_scale_bar,
        rendered_support_count=rendered_support_count,
        rendered_categorical_trait_count=rendered_categorical_trait_count,
        rendered_continuous_trait_count=0,
        rendered_metadata_strip_count=0,
        rendered_heatmap_column_count=0,
        collapsed_clade_count=0,
        missing_metadata_labels=sorted(set(missing_labels)),
    )
