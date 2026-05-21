from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class TreeRenderResult:
    """Stable report for one rendered phylogenetic tree SVG."""

    output_path: Path
    format: str
    layout: str
    tip_count: int
    visible_tip_count: int
    label_count: int
    has_scale_bar: bool
    scale_bar_length: float | None
    max_branch_distance: float | None
    rendered_support_count: int
    support_labels_validated: bool
    support_validation_warnings: list[str]
    rendered_categorical_trait_count: int
    rendered_continuous_trait_count: int
    rendered_metadata_strip_count: int
    rendered_heatmap_column_count: int
    rendered_internal_annotation_count: int
    rendered_branch_color_count: int
    rendered_internal_pie_count: int
    collapsed_clade_count: int
    missing_metadata_labels: list[str]


@dataclass(frozen=True, slots=True)
class AnnotationStrip:
    """Named taxon annotation column rendered beside a tree."""

    name: str
    values: dict[str, str]


@dataclass(frozen=True, slots=True)
class SupportLabelRenderAudit:
    """Decision about whether support labels are safe to render."""

    validated: bool
    labels_by_node: dict[str, str]
    warnings: list[str]
