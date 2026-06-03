from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.render.tree_svg import (
    SupportLabelRenderAudit,
    TreeRenderResult,
)


@dataclass(frozen=True, slots=True)
class FigureAnnotationCoverage:
    """Coverage and drift audit for one figure annotation surface."""

    surface: str
    aligned: bool
    covered_taxa: int
    missing_taxa: list[str]
    extra_taxa: list[str]


@dataclass(frozen=True, slots=True)
class FigureCollapsedCladeSummary:
    """Reviewer-facing summary for one collapsed clade."""

    clade_name: str
    descendant_count: int
    descendant_taxa: list[str]
    metadata_summaries: list[str]


@dataclass(frozen=True, slots=True)
class FigureLegendAudit:
    """Legend completeness check for rendered figure surfaces."""

    complete: bool
    entries: list[str]
    missing_entries: list[str]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class FigureLegendEntry:
    """One explicit legend entry emitted alongside a publication figure."""

    surface: str
    label: str
    swatch: str
    detail: str


@dataclass(frozen=True, slots=True)
class FigureLegibilityAudit:
    """Heuristic legibility review for one rendered publication tree figure."""

    legible: bool
    tip_label_font_size_px: int
    vertical_tip_spacing_px: int
    longest_visible_label_length: int
    estimated_longest_label_width_px: float
    available_label_lane_px: int
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class FigureCaptionDraft:
    """Structured caption draft for one publication-oriented tree figure."""

    title: str
    lead_sentence: str
    support_sentence: str
    scale_bar_sentence: str
    legend_sentence: str
    limitation_sentence: str
    caption_ready: bool


@dataclass(frozen=True, slots=True)
class FigureTableConsistencyReport:
    """Consistency check between rendered annotations and tabular exports."""

    consistent: bool
    missing_from_table: list[str]
    extra_in_table: list[str]
    label_mismatches: list[str]


@dataclass(slots=True)
class TreeFigureAuditReport:
    """Combined reviewer-facing audit for a rendered figure package."""

    support_audit: SupportLabelRenderAudit
    annotation_coverage: list[FigureAnnotationCoverage]
    collapsed_clades: list[FigureCollapsedCladeSummary]
    legend_audit: FigureLegendAudit
    table_consistency: FigureTableConsistencyReport
    scale_bar_valid: bool
    scale_bar_note: str
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class TreeFigurePackageResult:
    """Reviewer-facing artifact bundle for one rendered tree figure."""

    output_dir: Path
    figure_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    caption_path: Path
    legend_path: Path
    annotations_path: Path
    render: TreeRenderResult
    audit: TreeFigureAuditReport
    legend_entries: list[FigureLegendEntry]
    legibility_audit: FigureLegibilityAudit
    caption_draft: FigureCaptionDraft


__all__ = [
    "FigureAnnotationCoverage",
    "FigureCaptionDraft",
    "FigureCollapsedCladeSummary",
    "FigureLegendAudit",
    "FigureLegendEntry",
    "FigureLegibilityAudit",
    "FigureTableConsistencyReport",
    "TreeFigureAuditReport",
    "TreeFigurePackageResult",
]
