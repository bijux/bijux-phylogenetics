from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.biogeography.state_models import GeographicStateModelReport
from bijux_phylogenetics.phylogeography.geographic_map import GeographicMapReport
from bijux_phylogenetics.phylogeography.region_styles import (
    build_geographic_state_color_map,
    geographic_transition_support_colors,
    geographic_transition_support_details,
)
from bijux_phylogenetics.render.tree_svg import TreeRenderResult


@dataclass(frozen=True, slots=True)
class BiogeographyPublicationLegendEntry:
    """One explicit legend entry for publication-grade biogeography figures."""

    surface: str
    label: str
    swatch: str
    detail: str


@dataclass(frozen=True, slots=True)
class BiogeographyCaptionDraft:
    """Structured caption draft for a biogeography figure package."""

    title: str
    lead_sentence: str
    node_probability_sentence: str
    transition_sentence: str
    legend_sentence: str
    limitation_sentence: str
    caption_ready: bool


@dataclass(frozen=True, slots=True)
class BiogeographyPublicationAudit:
    """Reviewer-facing publication readiness audit for biogeography figures."""

    publication_ready: bool
    legend_complete: bool
    caption_ready: bool
    node_probabilities_visible: bool
    transitions_visible: bool
    map_state_colors_complete: bool
    rendered_internal_pie_count: int
    rendered_internal_probability_label_count: int
    expected_internal_node_count: int
    visible_transition_count: int
    state_color_count: int
    legend_entry_count: int
    reviewer_summary: list[str]
    limitations: list[str]


def build_biogeography_publication_legend_entries(
    *,
    state_report: GeographicStateModelReport,
    map_report: GeographicMapReport,
) -> list[BiogeographyPublicationLegendEntry]:
    """Build stable legend entries for tree and map publication surfaces."""
    state_labels = {
        row.most_likely_region
        for row in state_report.node_rows
        if row.most_likely_region
    }
    state_labels.update(
        row.state_label for row in map_report.marker_rows if row.state_label
    )
    state_colors = build_geographic_state_color_map(state_labels)
    transition_details = geographic_transition_support_details()
    transition_colors = geographic_transition_support_colors()
    entries = [
        BiogeographyPublicationLegendEntry(
            surface="state",
            label=label,
            swatch=color,
            detail="geographic region color is shared across the ancestral-state tree and map",
        )
        for label, color in state_colors.items()
    ]
    entries.extend(
        BiogeographyPublicationLegendEntry(
            surface="transition-support",
            label=label,
            swatch=transition_colors[label],
            detail=transition_details[label],
        )
        for label in ("strong", "moderate", "weak")
    )
    return entries


def write_biogeography_publication_legend(
    path: Path,
    entries: list[BiogeographyPublicationLegendEntry],
) -> Path:
    """Write the explicit legend ledger for a publication biogeography package."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["surface\tlabel\tswatch\tdetail"]
    lines.extend(
        "\t".join([entry.surface, entry.label, entry.swatch, entry.detail])
        for entry in entries
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def build_biogeography_caption_draft(
    *,
    state_report: GeographicStateModelReport,
    map_report: GeographicMapReport,
    audit: BiogeographyPublicationAudit | None = None,
) -> BiogeographyCaptionDraft:
    """Build a caption draft anchored in the rendered publication surfaces."""
    state_count = state_report.summary.observed_region_count
    caption_ready = audit.caption_ready if audit is not None else True
    limitations = (
        audit.limitations
        if audit is not None
        else ["publication review is incomplete until all figure audits pass"]
    )
    return BiogeographyCaptionDraft(
        title="Biogeographic ancestral-state tree and geographic transition map",
        lead_sentence=(
            f"Discrete ancestral-region reconstruction across {state_report.summary.analyzed_taxon_count} analyzed taxa resolves {state_count} observed regions on the rooted tree and linked geographic map."
        ),
        node_probability_sentence=(
            "Internal tree pies show node-wise region probabilities, and adjacent labels report the maximum posterior probability for each ancestral call."
        ),
        transition_sentence=(
            f"The map retains {map_report.summary.visible_line_count} visible geographic transitions, with line color keyed to strong, moderate, or weak support."
        ),
        legend_sentence=(
            "Geographic region colors are shared across the tree and map so state identity stays stable between ancestral nodes, tips, and geographic centroids."
        ),
        limitation_sentence=limitations[0],
        caption_ready=caption_ready,
    )


def write_biogeography_caption(
    path: Path,
    draft: BiogeographyCaptionDraft,
) -> Path:
    """Write one durable Markdown caption for biogeography publication figures."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"# {draft.title}",
                "",
                draft.lead_sentence,
                draft.node_probability_sentence,
                draft.transition_sentence,
                draft.legend_sentence,
                draft.limitation_sentence,
                "",
                f"caption_ready: {'true' if draft.caption_ready else 'false'}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return path


def build_biogeography_publication_audit(
    *,
    state_report: GeographicStateModelReport,
    map_report: GeographicMapReport,
    tree_render: TreeRenderResult,
    legend_entries: list[BiogeographyPublicationLegendEntry],
    exclusion_count: int,
) -> BiogeographyPublicationAudit:
    """Assess whether the biogeography figure package is publication-ready."""
    expected_internal_node_count = state_report.summary.internal_node_count
    node_probabilities_visible = (
        tree_render.rendered_internal_pie_count == expected_internal_node_count
        and tree_render.rendered_internal_annotation_count
        == expected_internal_node_count
    )
    state_labels = {
        row.most_likely_region
        for row in state_report.node_rows
        if row.most_likely_region
    }
    map_state_labels = {
        row.state_label for row in map_report.marker_rows if row.state_label
    }
    map_state_colors_complete = (
        state_labels.issubset(map_state_labels)
        and len(map_state_labels) >= state_report.summary.observed_region_count
    )
    transitions_visible = map_report.summary.visible_line_count > 0
    legend_complete = (
        len({entry.label for entry in legend_entries if entry.surface == "state"})
        >= state_report.summary.observed_region_count
    )
    caption_ready = node_probabilities_visible and map_state_colors_complete
    publication_ready = (
        node_probabilities_visible
        and transitions_visible
        and map_state_colors_complete
        and legend_complete
        and caption_ready
        and exclusion_count == 0
    )
    reviewer_summary = [
        f"internal probability pies rendered: {tree_render.rendered_internal_pie_count}/{expected_internal_node_count}",
        f"internal probability labels rendered: {tree_render.rendered_internal_annotation_count}/{expected_internal_node_count}",
        f"visible geographic transitions: {map_report.summary.visible_line_count}",
        f"geographic state colors represented on map: {len(map_state_labels)}/{state_report.summary.observed_region_count}",
    ]
    limitations: list[str] = []
    if not node_probabilities_visible:
        limitations.append(
            "internal node probabilities are not fully visible on the tree figure"
        )
    if not transitions_visible:
        limitations.append(
            "the geographic map does not show any visible region-transition lines"
        )
    if not map_state_colors_complete:
        limitations.append(
            "the geographic map does not represent every inferred region with a visible colored marker"
        )
    if exclusion_count:
        limitations.append(
            "one or more taxa, nodes, or transition events were excluded from the geographic figure surfaces"
        )
    if not limitations:
        limitations.append(
            "the current package is safe for publication-oriented geographic state review"
        )
    return BiogeographyPublicationAudit(
        publication_ready=publication_ready,
        legend_complete=legend_complete,
        caption_ready=caption_ready,
        node_probabilities_visible=node_probabilities_visible,
        transitions_visible=transitions_visible,
        map_state_colors_complete=map_state_colors_complete,
        rendered_internal_pie_count=tree_render.rendered_internal_pie_count,
        rendered_internal_probability_label_count=tree_render.rendered_internal_annotation_count,
        expected_internal_node_count=expected_internal_node_count,
        visible_transition_count=map_report.summary.visible_line_count,
        state_color_count=len(map_state_labels),
        legend_entry_count=len(legend_entries),
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )
