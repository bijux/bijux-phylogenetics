"""Rectangular and phylogram tree SVG layout rendering."""

from __future__ import annotations

from html import escape

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .contracts import AnnotationStrip
from .render_state import TreeSvgRenderState
from .shared import (
    Point,
    continuous_color,
    count_subtree_leaves,
    is_collapsed_node,
    node_signature,
    svg_pie_slices,
)
from .support_audit import coerce_support_label


def _node_x(
    *,
    layout: str,
    left_margin: int,
    horizontal_step: int,
    tree_width: float,
    max_distance: float,
    depth: int,
    distance: float,
) -> float:
    if layout == "phylogram" and max_distance > 0:
        return left_margin + (distance / max_distance) * tree_width
    return left_margin + depth * horizontal_step


def render_rectangular_tree_layout(
    tree: PhyloTree,
    *,
    layout: str,
    state: TreeSvgRenderState,
    labels: dict[str, str],
    collapsed_clade_names: set[str],
    row_height: int,
    top_margin: int,
    left_margin: int,
    horizontal_step: int,
    tree_width: float,
    max_distance: float,
    categorical_traits: dict[str, str],
    categorical_color_map: dict[str, str],
    continuous_traits: dict[str, float],
    continuous_min: float,
    continuous_max: float,
    metadata_strips: list[AnnotationStrip],
    metadata_strip_colors: list[dict[str, str]],
    heatmap_columns: list[AnnotationStrip],
    heatmap_specs: list[tuple[str, dict[str, str], float, float]],
    internal_annotations: dict[str, str],
    internal_annotation_colors: dict[str, str],
    branch_colors: dict[str, str],
    internal_pies: dict[str, dict[str, float]],
    internal_pie_colors: dict[str, str],
    show_support_values: bool,
    validated_support_labels: dict[str, str],
) -> None:
    """Render one rectangular or phylogram tree layout into the shared state."""

    def visit_rectangular(node: TreeNode, depth: int, distance: float) -> Point:
        branch_distance = distance + float(node.branch_length or 0.0)
        x = _node_x(
            layout=layout,
            left_margin=left_margin,
            horizontal_step=horizontal_step,
            tree_width=tree_width,
            max_distance=max_distance,
            depth=depth,
            distance=branch_distance if node is not tree.root else distance,
        )
        if node.is_leaf() or is_collapsed_node(node, collapsed_clade_names):
            y = top_margin + state.next_leaf_index * row_height + row_height / 2
            state.next_leaf_index += 1
            label = labels.get(node.name or "", node.name or "")
            if is_collapsed_node(node, collapsed_clade_names):
                label = (
                    f"{node.name or 'collapsed clade'}"
                    f" ({count_subtree_leaves(node)} tips)"
                )
                state.rendered_collapsed_clades += 1
                triangle = (
                    f"{x:.1f},{y - 11:.1f} {x + 24:.1f},{y:.1f} {x:.1f},{y + 11:.1f}"
                )
                state.overlays.append(
                    f'<polygon points="{triangle}" fill="#cbd5e1" stroke="#475569" stroke-width="1.5" class="collapsed-clade"/>'
                )
            elif node.name and node.name not in labels and labels:
                state.missing_labels.append(node.name)
            label_x = x + 18
            state.texts.append(
                f'<text x="{label_x:.1f}" y="{y + 5:.1f}" class="tip-label">{escape(label)}</text>'
            )
            trait_value = categorical_traits.get(node.name or "")
            if trait_value:
                marker_x = left_margin + tree_width + 190
                state.overlays.append(
                    f'<circle cx="{marker_x:.1f}" cy="{y - 4:.1f}" r="7" fill="{categorical_color_map[trait_value]}" class="trait-marker"/>'
                )
                state.rendered_categorical_trait_count += 1
            continuous_value = continuous_traits.get(node.name or "")
            if continuous_value is not None:
                bar_x = left_margin + tree_width + 210
                fill_fraction = (
                    1.0
                    if continuous_max <= continuous_min
                    else (continuous_value - continuous_min)
                    / (continuous_max - continuous_min)
                )
                fill_width = 52 * max(0.0, min(fill_fraction, 1.0))
                state.overlays.append(
                    f'<rect x="{bar_x:.1f}" y="{y - 11:.1f}" width="52" height="12" rx="6" fill="#e2e8f0" class="trait-bar-outline"/>'
                    f'<rect x="{bar_x:.1f}" y="{y - 11:.1f}" width="{fill_width:.1f}" height="12" rx="6" fill="{continuous_color(continuous_value, continuous_min, continuous_max)}" class="trait-bar-fill"/>'
                )
                state.rendered_continuous_trait_count += 1
            for strip_index, strip in enumerate(metadata_strips):
                strip_value = strip.values.get(node.name or "")
                if strip_value:
                    strip_x = left_margin + tree_width + 290 + strip_index * 24
                    state.overlays.append(
                        f'<rect x="{strip_x:.1f}" y="{y - 13:.1f}" width="16" height="16" rx="4" fill="{metadata_strip_colors[strip_index][strip_value]}" class="metadata-strip-cell"/>'
                    )
            heatmap_base_x = left_margin + tree_width + 318 + len(metadata_strips) * 24
            for column_index, column in enumerate(heatmap_columns):
                value = column.values.get(node.name or "")
                if not value:
                    continue
                cell_x = heatmap_base_x + column_index * 22
                mode, color_map, numeric_min, numeric_max = heatmap_specs[column_index]
                fill = (
                    continuous_color(float(value), numeric_min, numeric_max)
                    if mode == "numeric"
                    else color_map[value]
                )
                state.overlays.append(
                    f'<rect x="{cell_x:.1f}" y="{y - 13:.1f}" width="16" height="16" rx="4" fill="{fill}" class="heatmap-cell"/>'
                )
            return Point(x=x, y=y)

        child_points = [
            visit_rectangular(child, depth + 1, branch_distance)
            for child in node.children
        ]
        y = sum(point.y for point in child_points) / len(child_points)
        min_y = min(point.y for point in child_points)
        max_y = max(point.y for point in child_points)
        state.lines.append(
            f'<line x1="{x:.1f}" y1="{min_y:.1f}" x2="{x:.1f}" y2="{max_y:.1f}" class="branch"/>'
        )
        for child, child_point in zip(node.children, child_points, strict=False):
            branch_signature = node_signature(child)
            branch_color = branch_colors.get(branch_signature)
            branch_attributes = (
                f' stroke="{branch_color}" class="branch branch-colored"'
                if branch_color is not None
                else ' class="branch"'
            )
            state.lines.append(
                f'<line x1="{x:.1f}" y1="{child_point.y:.1f}" x2="{child_point.x:.1f}" y2="{child_point.y:.1f}"{branch_attributes}/>'
            )
            if branch_color is not None:
                state.rendered_branch_color_count += 1
        support_label = None
        if show_support_values:
            support_label = validated_support_labels.get(node_signature(node))
            if support_label is None:
                support_label = coerce_support_label(node.name)
        if support_label is not None and node is not tree.root:
            state.texts.append(
                f'<text x="{x + 8:.1f}" y="{y - 8:.1f}" class="support-label">{escape(support_label)}</text>'
            )
            state.rendered_support_count += 1
        annotation = internal_annotations.get(node_signature(node))
        pie_segments = internal_pies.get(node_signature(node))
        if pie_segments and not node.is_leaf():
            state.overlays.extend(
                svg_pie_slices(
                    x,
                    y,
                    7.0,
                    pie_segments,
                    internal_pie_colors,
                )
            )
            state.rendered_internal_pie_count += 1
        if annotation and not node.is_leaf():
            annotation_color = internal_annotation_colors.get(
                node_signature(node), "#7c3aed"
            )
            state.overlays.append(
                f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4.5" fill="{annotation_color}" class="internal-annotation-dot"/>'
            )
            state.texts.append(
                f'<text x="{x + 10:.1f}" y="{y + 4:.1f}" class="internal-annotation-label">{escape(annotation)}</text>'
            )
            state.rendered_internal_annotation_count += 1
        return Point(x=x, y=y)

    visit_rectangular(tree.root, 0, 0.0)
