from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.phylo.topology.tree import TreeNode


@dataclass(frozen=True, slots=True)
class Point:
    """One rendered SVG point in Cartesian coordinates."""

    x: float
    y: float


CATEGORICAL_PALETTE = (
    "#0f766e",
    "#1d4ed8",
    "#c2410c",
    "#7c3aed",
    "#b91c1c",
    "#047857",
    "#a16207",
    "#0f172a",
)


def count_subtree_leaves(node: TreeNode) -> int:
    """Count the total leaf descendants under one tree node."""
    if node.is_leaf():
        return 1
    return sum(count_subtree_leaves(child) for child in node.children)


def nice_scale_bar_length(max_distance: float) -> float:
    """Choose a stable reviewer-friendly scale-bar length."""
    if max_distance <= 0:
        return 0.0
    exponent = math.floor(math.log10(max_distance))
    base = 10**exponent
    for factor in (1.0, 0.5, 0.2, 0.1):
        candidate = base * factor
        if candidate <= max_distance / 3:
            return candidate
    return base / 10


def polar_point(
    center_x: float, center_y: float, radius: float, angle_radians: float
) -> Point:
    """Project one polar coordinate into SVG Cartesian coordinates."""
    return Point(
        x=center_x + radius * math.cos(angle_radians),
        y=center_y + radius * math.sin(angle_radians),
    )


def categorical_color_map(values: dict[str, str]) -> dict[str, str]:
    """Build a stable palette assignment for categorical annotation values."""
    categories = sorted({value for value in values.values() if value})
    return {
        category: CATEGORICAL_PALETTE[index % len(CATEGORICAL_PALETTE)]
        for index, category in enumerate(categories)
    }


def svg_pie_slices(
    center_x: float,
    center_y: float,
    radius: float,
    segments: dict[str, float],
    segment_colors: dict[str, str],
) -> list[str]:
    """Build deterministic SVG fragments for one internal node pie chart."""
    total = sum(value for value in segments.values() if value > 0.0)
    if total <= 0.0:
        return []
    start_angle = -math.pi / 2
    slices: list[str] = []
    for label, raw_value in segments.items():
        if raw_value <= 0.0:
            continue
        sweep = (raw_value / total) * 2 * math.pi
        end_angle = start_angle + sweep
        if sweep >= (2 * math.pi) - 1e-9:
            slices.append(
                f'<circle cx="{center_x:.1f}" cy="{center_y:.1f}" r="{radius:.1f}" fill="{segment_colors.get(label, "#94a3b8")}" class="internal-pie-slice"/>'
            )
        else:
            start_x = center_x + radius * math.cos(start_angle)
            start_y = center_y + radius * math.sin(start_angle)
            end_x = center_x + radius * math.cos(end_angle)
            end_y = center_y + radius * math.sin(end_angle)
            large_arc = 1 if sweep > math.pi else 0
            slices.append(
                f'<path d="M {center_x:.1f} {center_y:.1f} L {start_x:.1f} {start_y:.1f} A {radius:.1f} {radius:.1f} 0 {large_arc} 1 {end_x:.1f} {end_y:.1f} Z" fill="{segment_colors.get(label, "#94a3b8")}" class="internal-pie-slice"/>'
            )
        start_angle = end_angle
    return slices


def interpolate_channel(start: int, end: int, fraction: float) -> int:
    """Interpolate one RGB channel for continuous color ramps."""
    return round(start + (end - start) * fraction)


def continuous_color(value: float, minimum: float, maximum: float) -> str:
    """Map one continuous value into the house SVG gradient ramp."""
    fraction = 0.5 if maximum <= minimum else (value - minimum) / (maximum - minimum)
    red = interpolate_channel(219, 15, fraction)
    green = interpolate_channel(234, 118, fraction)
    blue = interpolate_channel(254, 110, fraction)
    return f"#{red:02x}{green:02x}{blue:02x}"


def is_numeric_strings(values: list[str]) -> bool:
    """Return whether a collection of strings is entirely numeric."""
    if not values:
        return False
    try:
        for value in values:
            float(value)
    except ValueError:
        return False
    return True


def is_collapsed_node(node: TreeNode, collapsed_clades: set[str]) -> bool:
    """Return whether a node should render as one collapsed clade marker."""
    return (
        not node.is_leaf() and node.name is not None and node.name in collapsed_clades
    )


def node_signature(node: TreeNode) -> str:
    """Build the stable render signature for one tree node."""
    if node.is_leaf():
        return node.name or "<unnamed>"
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(node_signature_taxa(child))
    return "|".join(sorted(taxa)) if taxa else node.name or "<unnamed>"


def node_signature_taxa(node: TreeNode) -> list[str]:
    """Collect descendant taxa used in one internal node signature."""
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(node_signature_taxa(child))
    return taxa


def count_visible_leaves(node: TreeNode, collapsed_clades: set[str]) -> int:
    """Count the leaves that remain visible after clade collapsing."""
    if node.is_leaf() or is_collapsed_node(node, collapsed_clades):
        return 1
    return sum(count_visible_leaves(child, collapsed_clades) for child in node.children)


def max_visible_depth(
    node: TreeNode, collapsed_clades: set[str], depth: int = 0
) -> int:
    """Compute the maximum visible topology depth after clade collapsing."""
    if node.is_leaf() or is_collapsed_node(node, collapsed_clades):
        return depth
    return max(
        max_visible_depth(child, collapsed_clades, depth + 1) for child in node.children
    )


def max_visible_distance(
    node: TreeNode, collapsed_clades: set[str], distance: float = 0.0
) -> float:
    """Compute the maximum visible branch-length distance after clade collapsing."""
    next_distance = distance + float(node.branch_length or 0.0)
    if node.is_leaf() or is_collapsed_node(node, collapsed_clades):
        return next_distance
    return max(
        max_visible_distance(child, collapsed_clades, next_distance)
        for child in node.children
    )
