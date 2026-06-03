"""Owned SVG rendering surface for phylogenetic trees."""

from __future__ import annotations

from .builder import render_tree_svg
from .contracts import AnnotationStrip, SupportLabelRenderAudit, TreeRenderResult
from .render_state import TreeSvgRenderState
from .shared import (
    Point,
    categorical_color_map,
    continuous_color,
    count_subtree_leaves,
    count_visible_leaves,
    is_collapsed_node,
    is_numeric_strings,
    max_visible_depth,
    max_visible_distance,
    nice_scale_bar_length,
    node_signature,
    node_signature_taxa,
    polar_point,
    svg_pie_slices,
)
from .support_audit import (
    audit_support_label_rendering,
    coerce_support_label,
    format_branch_value,
)

__all__ = [
    "AnnotationStrip",
    "Point",
    "SupportLabelRenderAudit",
    "TreeRenderResult",
    "TreeSvgRenderState",
    "audit_support_label_rendering",
    "categorical_color_map",
    "continuous_color",
    "coerce_support_label",
    "count_subtree_leaves",
    "count_visible_leaves",
    "format_branch_value",
    "is_collapsed_node",
    "is_numeric_strings",
    "max_visible_depth",
    "max_visible_distance",
    "nice_scale_bar_length",
    "node_signature",
    "node_signature_taxa",
    "polar_point",
    "render_tree_svg",
    "svg_pie_slices",
]
