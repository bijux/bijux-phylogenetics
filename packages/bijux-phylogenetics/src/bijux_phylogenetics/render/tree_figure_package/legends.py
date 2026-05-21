from __future__ import annotations

from bijux_phylogenetics.render.tree_svg import AnnotationStrip, TreeRenderResult

from .contracts import FigureLegendAudit, FigureLegendEntry


def build_legend_audit(
    *,
    render: TreeRenderResult,
    categorical_traits: dict[str, str],
    continuous_traits: dict[str, float],
    metadata_strips: list[AnnotationStrip],
    heatmap_columns: list[AnnotationStrip],
) -> FigureLegendAudit:
    """Check whether the figure legend fully documents rendered surfaces."""
    entries: list[str] = []
    if categorical_traits:
        entries.append("categorical trait legend")
    if continuous_traits:
        entries.append("continuous trait gradient")
    if metadata_strips:
        entries.extend(f"metadata strip: {strip.name}" for strip in metadata_strips)
    if heatmap_columns:
        entries.extend(f"heatmap column: {column.name}" for column in heatmap_columns)
    if render.rendered_support_count:
        entries.append("support label audit")
    if render.has_scale_bar:
        entries.append("branch-length scale bar")

    missing_entries: list[str] = []
    warnings: list[str] = []
    if render.rendered_support_count and not render.support_labels_validated:
        missing_entries.append("validated support label note")
        warnings.extend(render.support_validation_warnings)
    return FigureLegendAudit(
        complete=not missing_entries,
        entries=entries,
        missing_entries=missing_entries,
        warnings=warnings,
    )


def _categorical_legend_entries(
    *,
    surface: str,
    values: dict[str, str],
) -> list[FigureLegendEntry]:
    categories = sorted({value for value in values.values() if value})
    if not categories:
        return []
    colors = dict(
        zip(
            categories,
            (
                "#0f766e",
                "#1d4ed8",
                "#c2410c",
                "#7c3aed",
                "#b91c1c",
                "#047857",
                "#a16207",
                "#0f172a",
            ),
            strict=False,
        )
    )
    return [
        FigureLegendEntry(
            surface=surface,
            label=category,
            swatch=colors[category],
            detail=f"{surface} category rendered directly on the vector figure",
        )
        for category in categories
    ]


def _is_number(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def build_legend_entries(
    *,
    render: TreeRenderResult,
    categorical_traits: dict[str, str],
    continuous_traits: dict[str, float],
    metadata_strips: list[AnnotationStrip],
    heatmap_columns: list[AnnotationStrip],
) -> list[FigureLegendEntry]:
    """Build reviewer-facing explicit legend rows for one figure package."""
    entries: list[FigureLegendEntry] = []
    if render.has_scale_bar and render.scale_bar_length is not None:
        entries.append(
            FigureLegendEntry(
                surface="branch-length",
                label="scale bar",
                swatch="#0f172a",
                detail=f"scale bar represents branch length {render.scale_bar_length}",
            )
        )
    if render.rendered_support_count:
        entries.append(
            FigureLegendEntry(
                surface="support",
                label="validated support labels",
                swatch="#0f766e",
                detail=f"{render.rendered_support_count} support labels were rendered after audit",
            )
        )
    entries.extend(
        _categorical_legend_entries(
            surface="categorical trait",
            values=categorical_traits,
        )
    )
    if continuous_traits:
        minimum = min(continuous_traits.values())
        maximum = max(continuous_traits.values())
        entries.append(
            FigureLegendEntry(
                surface="continuous trait",
                label="gradient range",
                swatch=f"{minimum:.6g}..{maximum:.6g}",
                detail="continuous trait values are mapped to a low-to-high gradient",
            )
        )
    for strip in metadata_strips:
        entries.extend(
            _categorical_legend_entries(
                surface=f"metadata strip: {strip.name}",
                values=strip.values,
            )
        )
    for column in heatmap_columns:
        observed_values = [value for value in column.values.values() if value]
        if not observed_values:
            continue
        if all(_is_number(value) for value in observed_values):
            numeric_values = [float(value) for value in observed_values]
            entries.append(
                FigureLegendEntry(
                    surface=f"heatmap column: {column.name}",
                    label="numeric gradient",
                    swatch=f"{min(numeric_values):.6g}..{max(numeric_values):.6g}",
                    detail="heatmap column uses a continuous low-to-high color gradient",
                )
            )
        else:
            entries.extend(
                _categorical_legend_entries(
                    surface=f"heatmap column: {column.name}",
                    values=column.values,
                )
            )
    return entries


__all__ = [
    "build_legend_audit",
    "build_legend_entries",
]
