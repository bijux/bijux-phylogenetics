from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ...continuous import BrownianMotionFitReport, OUTraitModelReport
from ...model_selection import ComparativeModelComparisonReport


@dataclass(frozen=True, slots=True)
class ComparativeModelFigureLegendEntry:
    """One legend row for the comparative model-comparison package."""

    surface: str
    label: str
    swatch: str
    detail: str


@dataclass(frozen=True, slots=True)
class ComparativeModelFigureCaptionDraft:
    """Structured caption draft for one comparative model-comparison package."""

    title: str
    lead_sentence: str
    criteria_sentence: str
    likelihood_sentence: str
    parameter_sentence: str
    fit_sentence: str
    limitation_sentence: str
    caption_ready: bool


@dataclass(frozen=True, slots=True)
class ComparativeModelCriteriaRow:
    """Reviewer-facing information-criterion summary for one fitted model."""

    model: str
    parameter_count: int
    log_likelihood: float
    aic: float
    aicc: float
    delta_aic: float
    delta_aicc: float
    selected: bool


@dataclass(frozen=True, slots=True)
class ComparativeModelLikelihoodRow:
    """Reviewer-facing likelihood summary for one fitted model."""

    model: str
    log_likelihood: float
    delta_log_likelihood: float
    selected: bool


@dataclass(frozen=True, slots=True)
class ComparativeModelParameterRow:
    """Reviewer-facing parameter summary for one fitted model parameter."""

    model: str
    parameter: str
    estimate: float
    lower_95: float
    upper_95: float
    interval_method: str


@dataclass(frozen=True, slots=True)
class ComparativeModelFitRow:
    """Reviewer-facing fit diagnostics for one fitted comparative model."""

    model: str
    taxon_count: int
    residual_variance: float
    max_abs_standardized_residual: float
    phylogenetic_residual_lambda: float
    outlier_taxon_count: int
    warning_count: int
    convergence_status: str
    selected: bool


@dataclass(frozen=True, slots=True)
class ComparativeModelFigureAudit:
    """Publication-oriented audit for a comparative model-comparison package."""

    publication_ready: bool
    criteria_surface_visible: bool
    likelihood_surface_visible: bool
    parameter_surface_visible: bool
    fit_surface_visible: bool
    legend_complete: bool
    caption_ready: bool
    finite_aicc_model_count: int
    support_distinct: bool
    selected_model: str
    aicc_delta: float | None
    plotted_model_count: int
    rendered_parameter_count: int
    rendered_fit_row_count: int
    warning_count: int
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class ComparativeModelFigurePackageResult:
    output_dir: Path
    criteria_figure_path: Path
    likelihood_figure_path: Path
    parameter_figure_path: Path
    fit_figure_path: Path
    criteria_table_path: Path
    likelihood_table_path: Path
    parameter_table_path: Path
    fit_table_path: Path
    legend_path: Path
    caption_path: Path
    review_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    comparison_report: ComparativeModelComparisonReport
    brownian_report: BrownianMotionFitReport
    ou_report: OUTraitModelReport
    criteria_rows: list[ComparativeModelCriteriaRow]
    likelihood_rows: list[ComparativeModelLikelihoodRow]
    parameter_rows: list[ComparativeModelParameterRow]
    fit_rows: list[ComparativeModelFitRow]
    legend_entries: list[ComparativeModelFigureLegendEntry]
    caption_draft: ComparativeModelFigureCaptionDraft
    audit: ComparativeModelFigureAudit
    machine_manifest: dict[str, object]
