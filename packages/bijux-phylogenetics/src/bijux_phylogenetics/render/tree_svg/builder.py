"""Build deterministic SVG tree figures from one phylogenetic tree input."""

from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.diagnostics.validation.structure import _load_tree

from .circular_layout import render_circular_tree_layout
from .contracts import AnnotationStrip, TreeRenderResult
from .presentation import build_tree_svg_document
from .rectangular_layout import render_rectangular_tree_layout
from .render_state import TreeSvgRenderState
from .shared import (
    categorical_color_map as build_categorical_color_map,
)
from .shared import (
    count_visible_leaves,
    is_numeric_strings,
    max_visible_depth,
    max_visible_distance,
)


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

    if layout == "circular":
        render_circular_tree_layout(
            tree,
            state=state,
            width=width,
            height=height,
            max_depth=max_depth,
            max_distance=max_distance,
            leaf_count=leaf_count,
            labels=labels,
            collapsed_clade_names=collapsed_clade_names,
            categorical_traits=categorical_traits,
            categorical_color_map=categorical_color_map,
            continuous_traits=continuous_traits,
            continuous_min=continuous_min,
            continuous_max=continuous_max,
            internal_annotations=internal_annotations,
            internal_annotation_colors=internal_annotation_colors,
            branch_colors=branch_colors,
            internal_pies=internal_pies,
            internal_pie_colors=internal_pie_colors,
            show_support_values=show_support_values,
            validated_support_labels=validated_support_labels,
        )
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

    svg, has_scale_bar, scale_bar_length = build_tree_svg_document(
        width=width,
        height=height,
        layout=layout,
        max_distance=max_distance,
        tree_width=tree_width,
        left_margin=left_margin,
        top_margin=top_margin,
        metadata_strips=metadata_strips,
        heatmap_columns=heatmap_columns,
        categorical_color_map=categorical_color_map,
        continuous_values=continuous_values,
        continuous_min=continuous_min,
        continuous_max=continuous_max,
        state=state,
    )
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
