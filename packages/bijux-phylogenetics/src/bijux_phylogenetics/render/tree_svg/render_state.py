"""Mutable rendering state owned by one tree SVG build."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(slots=True)
class TreeSvgRenderState:
    """Accumulate SVG fragments and render counters for one tree image."""

    lines: list[str] = field(default_factory=list)
    texts: list[str] = field(default_factory=list)
    overlays: list[str] = field(default_factory=list)
    missing_labels: list[str] = field(default_factory=list)
    next_leaf_index: int = 0
    rendered_support_count: int = 0
    rendered_categorical_trait_count: int = 0
    rendered_continuous_trait_count: int = 0
    rendered_internal_annotation_count: int = 0
    rendered_branch_color_count: int = 0
    rendered_internal_pie_count: int = 0
    rendered_collapsed_clades: int = 0
