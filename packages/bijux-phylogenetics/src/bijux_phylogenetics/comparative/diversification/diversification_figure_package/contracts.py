from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..models import (
    CladeDiversificationScanReport,
    DiversificationMethodReport,
    DiversificationMethodsSummaryTextResult,
    DiversificationModelComparisonReport,
    LineageThroughTimeReport,
    SamplingFractionReport,
)


@dataclass(frozen=True, slots=True)
class DiversificationFigureLegendEntry:
    """One legend row for the diversification figure package."""

    surface: str
    label: str
    swatch: str
    detail: str


@dataclass(frozen=True, slots=True)
class DiversificationFigureCaptionDraft:
    """Structured caption draft for one diversification figure package."""

    title: str
    lead_sentence: str
    lineage_sentence: str
    clade_sentence: str
    model_sentence: str
    limitation_sentence: str
    caption_ready: bool


@dataclass(frozen=True, slots=True)
class DiversificationFigureAudit:
    """Publication-oriented audit for diversification figure completeness."""

    publication_ready: bool
    lineage_curve_visible: bool
    clade_outlier_surface_visible: bool
    model_comparison_visible: bool
    legend_complete: bool
    caption_ready: bool
    sampling_metadata_complete: bool | None
    plotted_ltt_point_count: int
    plotted_clade_count: int
    highlighted_outlier_count: int
    plotted_model_count: int
    better_model: str
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class DiversificationFigurePackageResult:
    output_dir: Path
    lineage_figure_path: Path
    clade_figure_path: Path
    model_figure_path: Path
    lineage_table_path: Path
    clade_table_path: Path
    model_table_path: Path
    legend_path: Path
    caption_path: Path
    methods_summary_path: Path
    review_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    lineage_report: LineageThroughTimeReport
    clade_report: CladeDiversificationScanReport
    model_report: DiversificationModelComparisonReport
    sampling_report: SamplingFractionReport | None
    methods_report: DiversificationMethodReport
    methods_summary: DiversificationMethodsSummaryTextResult
    legend_entries: list[DiversificationFigureLegendEntry]
    caption_draft: DiversificationFigureCaptionDraft
    audit: DiversificationFigureAudit
    machine_manifest: dict[str, object]
