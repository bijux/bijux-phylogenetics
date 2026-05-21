"""Compatibility facade for the owned tree SVG rendering package."""

from __future__ import annotations

from bijux_phylogenetics.render.tree_svg import (
    AnnotationStrip,
    SupportLabelRenderAudit,
    TreeRenderResult,
    audit_support_label_rendering,
    render_tree_svg,
)

__all__ = [
    "AnnotationStrip",
    "SupportLabelRenderAudit",
    "TreeRenderResult",
    "audit_support_label_rendering",
    "render_tree_svg",
]
