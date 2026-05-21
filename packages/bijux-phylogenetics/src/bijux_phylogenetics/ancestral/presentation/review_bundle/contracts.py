from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.continuous import (
    ContinuousAncestralExclusion,
    ContinuousAncestralReport,
    ContinuousAncestralSummary,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralExclusion,
    DiscreteAncestralReport,
    DiscreteAncestralSummary,
)
from bijux_phylogenetics.ancestral.discrete.review import AncestralTransitionReport
from bijux_phylogenetics.ancestral.presentation.methods_text import (
    AncestralMethodsSummaryTextResult,
)
from bijux_phylogenetics.ancestral.presentation.visualization import (
    AncestralVisualizationResult,
)
from bijux_phylogenetics.reports.review import ReviewerAuditChecklist


@dataclass(frozen=True, slots=True)
class AncestralContinuousChangeBranchRow:
    parent_node: str
    child_node: str
    child_descendant_taxa: tuple[str, ...]
    branch_length: float | None
    parent_estimate: float
    child_estimate: float
    delta: float
    absolute_delta: float
    direction: str


@dataclass(frozen=True, slots=True)
class AncestralContinuousChangeCountRow:
    direction: str
    branch_count: int
    branch_fraction: float
    mean_delta: float
    minimum_delta: float
    maximum_delta: float


@dataclass(slots=True)
class AncestralReportPackageResult:
    output_dir: Path
    report_path: Path
    methods_summary_path: Path
    reviewer_audit_checklist_path: Path
    figure_path: Path
    figure_png_path: Path
    figure_html_path: Path
    summary_table_path: Path
    node_table_path: Path
    uncertainty_table_path: Path
    transition_count_table_path: Path
    transition_branch_table_path: Path
    exclusion_table_path: Path
    manifest_path: Path
    reconstruction_kind: str
    model: str
    methods_summary: AncestralMethodsSummaryTextResult
    summary: ContinuousAncestralSummary | DiscreteAncestralSummary
    reconstruction: ContinuousAncestralReport | DiscreteAncestralReport
    figure: AncestralVisualizationResult
    figure_png: AncestralVisualizationResult
    figure_html: AncestralVisualizationResult
    transition_count_rows: list[AncestralContinuousChangeCountRow]
    transition_branch_rows: list[AncestralContinuousChangeBranchRow]
    transition_report: AncestralTransitionReport | None
    exclusions: list[ContinuousAncestralExclusion | DiscreteAncestralExclusion]
    reviewer_audit_checklist: ReviewerAuditChecklist
    machine_manifest: dict[str, object]
