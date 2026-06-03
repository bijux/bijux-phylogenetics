"""Circular tree SVG layout rendering."""

from __future__ import annotations

from html import escape
import math

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .render_state import TreeSvgRenderState
from .shared import (
    continuous_color,
    count_subtree_leaves,
    is_collapsed_node,
    node_signature,
    polar_point,
    svg_pie_slices,
)
from .support_audit import coerce_support_label


def render_circular_tree_layout(
    tree: PhyloTree,
    *,
    state: TreeSvgRenderState,
    width: int,
    height: int,
    max_depth: int,
    max_distance: float,
    leaf_count: int,
    labels: dict[str, str],
    collapsed_clade_names: set[str],
    categorical_traits: dict[str, str],
    categorical_color_map: dict[str, str],
    continuous_traits: dict[str, float],
    continuous_min: float,
    continuous_max: float,
    internal_annotations: dict[str, str],
    internal_annotation_colors: dict[str, str],
    branch_colors: dict[str, str],
    internal_pies: dict[str, dict[str, float]],
    internal_pie_colors: dict[str, str],
    show_support_values: bool,
    validated_support_labels: dict[str, str],
) -> None:
    """Render one circular tree layout into the shared SVG state."""

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
            angle = (2 * math.pi * state.next_leaf_index / max(leaf_count, 1)) - (
                math.pi / 2
            )
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
                left_point = polar_point(center_x, center_y, radial + 18, angle - 0.08)
                right_point = polar_point(center_x, center_y, radial + 18, angle + 0.08)
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
            child_branch_distance = branch_distance + float(child.branch_length or 0.0)
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
            support_point = polar_point(center_x, center_y, radial + 14, node_angle)
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
            annotation_point = polar_point(center_x, center_y, radial + 12, node_angle)
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
