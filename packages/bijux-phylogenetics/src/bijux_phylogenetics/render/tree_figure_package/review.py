from __future__ import annotations

from bijux_phylogenetics.render.tree_svg import TreeRenderResult

from .contracts import (
    FigureCaptionDraft,
    FigureLegendEntry,
    FigureLegibilityAudit,
    TreeFigureAuditReport,
)


def build_legibility_audit(
    *,
    render: TreeRenderResult,
    visible_labels: list[str],
    has_annotation_columns: bool,
) -> FigureLegibilityAudit:
    """Evaluate whether the rendered figure still fits a publication lane."""
    tip_label_font_size_px = 16
    vertical_tip_spacing_px = 56 if render.layout != "circular" else 20
    longest_visible_label_length = max(
        (len(label) for label in visible_labels), default=0
    )
    estimated_longest_label_width_px = round(
        longest_visible_label_length * tip_label_font_size_px * 0.58,
        2,
    )
    available_label_lane_px = (
        172 if has_annotation_columns and render.layout != "circular" else 420
    )
    warnings: list[str] = []
    if estimated_longest_label_width_px > available_label_lane_px:
        warnings.append(
            "the longest visible label likely exceeds the reserved label lane and should be shortened or moved to metadata"
        )
    if (
        render.layout != "circular"
        and vertical_tip_spacing_px < tip_label_font_size_px * 2
    ):
        warnings.append(
            "vertical tip spacing fell below the publication legibility threshold"
        )
    if render.visible_tip_count > 80:
        warnings.append(
            "the visible tip count is high enough that clade collapsing or panel subdivision should be considered before journal submission"
        )
    return FigureLegibilityAudit(
        legible=not warnings,
        tip_label_font_size_px=tip_label_font_size_px,
        vertical_tip_spacing_px=vertical_tip_spacing_px,
        longest_visible_label_length=longest_visible_label_length,
        estimated_longest_label_width_px=estimated_longest_label_width_px,
        available_label_lane_px=available_label_lane_px,
        warnings=warnings,
    )


def build_caption_draft(
    *,
    title: str,
    render: TreeRenderResult,
    audit: TreeFigureAuditReport,
    legend_entries: list[FigureLegendEntry],
) -> FigureCaptionDraft:
    """Draft a publication caption from the rendered figure and reviewer audit."""
    lead_sentence = f"{title} shows {render.visible_tip_count} rendered taxa from a source tree with {render.tip_count} total tips using a {render.layout} layout."
    support_sentence = (
        f"Validated branch support labels are shown for {render.rendered_support_count} internal branches."
        if render.rendered_support_count
        else "Branch support labels were omitted because no validated support surface was available."
    )
    scale_bar_sentence = (
        f"Branch lengths are scaled directly on the figure using a scale bar of {render.scale_bar_length}."
        if render.has_scale_bar and render.scale_bar_length is not None
        else "A branch-length scale bar is not shown because the selected layout is not branch-length proportional."
    )
    legend_sentence = f"The figure legend contains {len(legend_entries)} explicit entries covering rendered trait, metadata, support, and scale surfaces."
    limitation_sentence = (
        "Reviewer-facing audits did not record additional publication limitations."
        if not (audit.limitations or audit.legend_audit.warnings)
        else "Publication review still requires attention to: "
        + "; ".join([*audit.limitations, *audit.legend_audit.warnings])
        + "."
    )
    return FigureCaptionDraft(
        title=title,
        lead_sentence=lead_sentence,
        support_sentence=support_sentence,
        scale_bar_sentence=scale_bar_sentence,
        legend_sentence=legend_sentence,
        limitation_sentence=limitation_sentence,
        caption_ready=audit.scale_bar_valid
        and audit.legend_audit.complete
        and audit.table_consistency.consistent,
    )


__all__ = [
    "build_caption_draft",
    "build_legibility_audit",
]
