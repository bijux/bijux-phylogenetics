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


def _count_subtree_leaves(node: TreeNode) -> int:
    if node.is_leaf():
        return 1
    return sum(_count_subtree_leaves(child) for child in node.children)


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


def _interpolate_channel(start: int, end: int, fraction: float) -> int:
    return round(start + (end - start) * fraction)


def _continuous_color(value: float, minimum: float, maximum: float) -> str:
    if maximum <= minimum:
        fraction = 0.5
    else:
        fraction = (value - minimum) / (maximum - minimum)
    red = _interpolate_channel(219, 15, fraction)
    green = _interpolate_channel(234, 118, fraction)
    blue = _interpolate_channel(254, 110, fraction)
    return f"#{red:02x}{green:02x}{blue:02x}"


def _is_collapsed_node(node: TreeNode, collapsed_clades: set[str]) -> bool:
    return not node.is_leaf() and node.name is not None and node.name in collapsed_clades


def _count_visible_leaves(node: TreeNode, collapsed_clades: set[str]) -> int:
    if node.is_leaf() or _is_collapsed_node(node, collapsed_clades):
        return 1
    return sum(_count_visible_leaves(child, collapsed_clades) for child in node.children)


def _max_visible_depth(node: TreeNode, collapsed_clades: set[str], depth: int = 0) -> int:
    if node.is_leaf() or _is_collapsed_node(node, collapsed_clades):
        return depth
    return max(_max_visible_depth(child, collapsed_clades, depth + 1) for child in node.children)


def _max_visible_distance(node: TreeNode, collapsed_clades: set[str], distance: float = 0.0) -> float:
    next_distance = distance + float(node.branch_length or 0.0)
    if node.is_leaf() or _is_collapsed_node(node, collapsed_clades):
        return next_distance
    return max(_max_visible_distance(child, collapsed_clades, next_distance) for child in node.children)


def render_tree_svg(
    tree_path: Path,
    *,
    out_path: Path,
    labels: dict[str, str] | None = None,
    layout: str = "cladogram",
    show_support_values: bool = False,
    categorical_traits: dict[str, str] | None = None,
    continuous_traits: dict[str, float] | None = None,
    collapsed_clades: list[str] | None = None,
) -> TreeRenderResult:
    """Render a deterministic SVG tree as a cladogram, phylogram, or circular tree."""
    if layout not in {"cladogram", "phylogram", "circular"}:
        raise ValueError(f"unsupported tree layout: {layout}")

    tree = _load_tree(tree_path)
    labels = labels or {}
    collapsed_clade_names = {name for name in (collapsed_clades or []) if name}
    row_height = 56
    left_margin = 48
    right_margin = 320
    top_margin = 40
    bottom_margin = 72
    horizontal_step = 150
    scale_width = 520
    leaf_count = _count_visible_leaves(tree.root, collapsed_clade_names)
    max_depth = max(_max_visible_depth(tree.root, collapsed_clade_names), 1)
    max_distance = (
        _max_visible_distance(tree.root, collapsed_clade_names, 0.0)
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

    lines: list[str] = []
    texts: list[str] = []
    overlays: list[str] = []
    missing_labels: list[str] = []
    next_leaf_index = 0
    rendered_support_count = 0
    rendered_categorical_trait_count = 0
    rendered_continuous_trait_count = 0
    rendered_collapsed_clades = 0
    categorical_traits = categorical_traits or {}
    categorical_color_map = _categorical_color_map(categorical_traits)
    continuous_traits = continuous_traits or {}
    continuous_values = list(continuous_traits.values())
    continuous_min = min(continuous_values) if continuous_values else 0.0
    continuous_max = max(continuous_values) if continuous_values else 0.0

    def node_x(depth: int, distance: float) -> float:
        if layout == "phylogram" and max_distance > 0:
            return left_margin + (distance / max_distance) * tree_width
        return left_margin + depth * horizontal_step

    def visit_rectangular(node: TreeNode, depth: int, distance: float) -> _Point:
        nonlocal next_leaf_index
        nonlocal rendered_support_count
        nonlocal rendered_categorical_trait_count
        nonlocal rendered_continuous_trait_count
        nonlocal rendered_collapsed_clades
        branch_distance = distance + float(node.branch_length or 0.0)
        x = node_x(depth, branch_distance if node is not tree.root else distance)
        if node.is_leaf() or _is_collapsed_node(node, collapsed_clade_names):
            y = top_margin + next_leaf_index * row_height + row_height / 2
            next_leaf_index += 1
            label = labels.get(node.name or "", node.name or "")
            if _is_collapsed_node(node, collapsed_clade_names):
                label = f"{node.name or 'collapsed clade'} ({_count_subtree_leaves(node)} tips)"
                rendered_collapsed_clades += 1
                triangle = (
                    f"{x:.1f},{y - 11:.1f} "
                    f"{x + 24:.1f},{y:.1f} "
                    f"{x:.1f},{y + 11:.1f}"
                )
                overlays.append(
                    f'<polygon points="{triangle}" fill="#cbd5e1" stroke="#475569" stroke-width="1.5" class="collapsed-clade"/>'
                )
            elif node.name and node.name not in labels and labels:
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
            continuous_value = continuous_traits.get(node.name or "")
            if continuous_value is not None:
                bar_x = left_margin + tree_width + 210
                fill_fraction = 1.0 if continuous_max <= continuous_min else (continuous_value - continuous_min) / (continuous_max - continuous_min)
                fill_width = 52 * max(0.0, min(fill_fraction, 1.0))
                overlays.append(
                    f'<rect x="{bar_x:.1f}" y="{y - 11:.1f}" width="52" height="12" rx="6" fill="#e2e8f0" class="trait-bar-outline"/>'
                    f'<rect x="{bar_x:.1f}" y="{y - 11:.1f}" width="{fill_width:.1f}" height="12" rx="6" fill="{_continuous_color(continuous_value, continuous_min, continuous_max)}" class="trait-bar-fill"/>'
                )
                rendered_continuous_trait_count += 1
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
            if node.is_leaf() or _is_collapsed_node(node, collapsed_clade_names):
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
            nonlocal rendered_continuous_trait_count
            nonlocal rendered_collapsed_clades
            branch_distance = distance + float(node.branch_length or 0.0)
            radial = radial_distance(depth, branch_distance if node is not tree.root else distance)
            if node.is_leaf() or _is_collapsed_node(node, collapsed_clade_names):
                angle = angle_cache[id(node)]
                label = labels.get(node.name or "", node.name or "")
                if _is_collapsed_node(node, collapsed_clade_names):
                    label = f"{node.name or 'collapsed clade'} ({_count_subtree_leaves(node)} tips)"
                    rendered_collapsed_clades += 1
                elif node.name and node.name not in labels and labels:
                    missing_labels.append(node.name)
                anchor = "start" if math.cos(angle) >= 0 else "end"
                label_point = _polar_point(center_x, center_y, radial + 18, angle)
                texts.append(
                    f'<text x="{label_point.x:.1f}" y="{label_point.y + 5:.1f}" text-anchor="{anchor}" class="tip-label">{escape(label)}</text>'
                )
                if _is_collapsed_node(node, collapsed_clade_names):
                    node_point = _polar_point(center_x, center_y, radial, angle)
                    left_point = _polar_point(center_x, center_y, radial + 18, angle - 0.08)
                    right_point = _polar_point(center_x, center_y, radial + 18, angle + 0.08)
                    overlays.append(
                        f'<polygon points="{node_point.x:.1f},{node_point.y:.1f} {left_point.x:.1f},{left_point.y:.1f} {right_point.x:.1f},{right_point.y:.1f}" fill="#cbd5e1" stroke="#475569" stroke-width="1.5" class="collapsed-clade"/>'
                    )
                trait_value = categorical_traits.get(node.name or "")
                if trait_value:
                    marker_point = _polar_point(center_x, center_y, radial + 10, angle)
                    overlays.append(
                        f'<circle cx="{marker_point.x:.1f}" cy="{marker_point.y:.1f}" r="6" fill="{categorical_color_map[trait_value]}" class="trait-marker"/>'
                    )
                    rendered_categorical_trait_count += 1
                continuous_value = continuous_traits.get(node.name or "")
                if continuous_value is not None:
                    marker_point = _polar_point(center_x, center_y, radial + 10, angle)
                    overlays.append(
                        f'<circle cx="{marker_point.x:.1f}" cy="{marker_point.y:.1f}" r="6" fill="{_continuous_color(continuous_value, continuous_min, continuous_max)}" class="trait-marker"/>'
                    )
                    rendered_continuous_trait_count += 1
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

    continuous_legend = ""
    if continuous_values:
        legend_x = width - 220
        legend_y = height - 60
        continuous_legend = (
            f'<text x="{legend_x - 2:.1f}" y="{legend_y - 10:.1f}" class="legend-title">continuous trait</text>'
            f'<defs><linearGradient id="continuous-trait-gradient" x1="0%" y1="0%" x2="100%" y2="0%">'
            f'<stop offset="0%" stop-color="{_continuous_color(continuous_min, continuous_min, continuous_max)}"/>'
            f'<stop offset="100%" stop-color="{_continuous_color(continuous_max, continuous_min, continuous_max)}"/>'
            f'</linearGradient></defs>'
            f'<rect x="{legend_x:.1f}" y="{legend_y:.1f}" width="120" height="12" rx="6" fill="url(#continuous-trait-gradient)"/>'
            f'<text x="{legend_x:.1f}" y="{legend_y + 28:.1f}" class="legend-label">{escape(_format_branch_value(continuous_min))}</text>'
            f'<text x="{legend_x + 88:.1f}" y="{legend_y + 28:.1f}" class="legend-label">{escape(_format_branch_value(continuous_max))}</text>'
        )

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="phylogenetic tree">
  <style>
    .panel {{ fill: #f7fbfa; stroke: #d7e3e1; stroke-width: 1; rx: 18; ry: 18; }}
    .branch {{ stroke: #0f172a; stroke-width: 2.2; stroke-linecap: round; fill: none; }}
    .scale-bar {{ stroke: #0f172a; stroke-width: 2; stroke-linecap: round; }}
    .collapsed-clade {{ fill-opacity: 0.95; }}
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
  {continuous_legend}
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
        rendered_continuous_trait_count=rendered_continuous_trait_count,
        rendered_metadata_strip_count=0,
        rendered_heatmap_column_count=0,
        collapsed_clade_count=rendered_collapsed_clades,
        missing_metadata_labels=sorted(set(missing_labels)),
    )
