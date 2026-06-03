from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    AlignmentForensicReport,
    AlignmentSummary,
    AlignmentWindowSummary,
)
from bijux_phylogenetics.reports.review import ReviewerAuditChecklist


@dataclass(frozen=True, slots=True)
class AlignmentHeatmapCell:
    """One sequence-by-site-bin uncertainty cell for the alignment heatmap."""

    identifier: str
    bin_start: int
    bin_end: int
    uncertainty_fraction: float
    gap_fraction: float
    missing_fraction: float
    ambiguity_fraction: float


@dataclass(frozen=True, slots=True)
class AlignmentFigureLegendEntry:
    """One explicit legend entry for the alignment figure package."""

    surface: str
    label: str
    swatch: str
    detail: str


@dataclass(frozen=True, slots=True)
class AlignmentFigureCaptionDraft:
    """Structured caption draft for the alignment figure package."""

    title: str
    lead_sentence: str
    heatmap_sentence: str
    site_summary_sentence: str
    sequence_panel_sentence: str
    limitation_sentence: str
    caption_ready: bool


@dataclass(frozen=True, slots=True)
class AlignmentFigureAudit:
    """Publication-oriented audit for alignment figure completeness and risk."""

    publication_ready: bool
    heatmap_visible: bool
    site_summary_visible: bool
    sequence_panel_visible: bool
    legend_complete: bool
    caption_ready: bool
    suspicious_alignment: bool
    quality_score: float
    heatmap_row_count: int
    heatmap_bin_count: int
    plotted_window_count: int
    plotted_sequence_count: int
    invalid_character_count: int
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class AlignmentFigurePackageResult:
    output_dir: Path
    heatmap_figure_path: Path
    site_summary_figure_path: Path
    sequence_panel_figure_path: Path
    heatmap_table_path: Path
    window_table_path: Path
    ranking_table_path: Path
    legend_path: Path
    caption_path: Path
    review_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    reviewer_audit_checklist_path: Path
    summary: AlignmentSummary
    forensic: AlignmentForensicReport
    windows: list[AlignmentWindowSummary]
    heatmap_cells: list[AlignmentHeatmapCell]
    legend_entries: list[AlignmentFigureLegendEntry]
    caption_draft: AlignmentFigureCaptionDraft
    audit: AlignmentFigureAudit
    reviewer_audit_checklist: ReviewerAuditChecklist
    machine_manifest: dict[str, object]
