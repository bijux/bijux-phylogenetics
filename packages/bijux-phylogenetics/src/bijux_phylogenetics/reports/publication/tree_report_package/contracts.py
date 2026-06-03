from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import (
    TreeForensicReport,
    TreeInspectionReport,
    TreeValidationReport,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import MethodTierAssessment
from bijux_phylogenetics.render.tree_svg import (
    SupportLabelRenderAudit,
    TreeRenderResult,
)
from bijux_phylogenetics.reports.methods import (
    TreeValidationMethodsSummaryTextResult,
)
from bijux_phylogenetics.reports.review import ReviewerAuditChecklist
from bijux_phylogenetics.trees import (
    BranchLengthDistributionReport,
    CladeTableReport,
)


@dataclass(frozen=True, slots=True)
class TreeSupportRow:
    node_kind: str
    node: str
    node_label: str | None
    descendant_taxa: tuple[str, ...]
    support: float | None
    support_fraction: float | None
    support_class: str
    branch_length: float | None
    root_depth: float | None


@dataclass(frozen=True, slots=True)
class TreeBranchStatisticsRow:
    branch_count: int
    defined_branch_count: int
    missing_branch_count: int
    zero_length_branch_count: int
    negative_branch_count: int
    positive_branch_count: int
    long_outlier_count: int
    short_outlier_count: int
    minimum_branch_length: float | None
    maximum_branch_length: float | None
    mean_branch_length: float | None
    median_branch_length: float | None
    positive_branch_median: float | None


@dataclass(slots=True)
class TreeReportPackageResult:
    output_dir: Path
    report_path: Path
    figure_path: Path
    methods_summary_path: Path
    reviewer_audit_checklist_path: Path
    support_table_path: Path
    clade_table_path: Path
    branch_stats_path: Path
    manifest_path: Path
    validation: TreeValidationReport
    inspection: TreeInspectionReport
    forensic: TreeForensicReport
    figure: TreeRenderResult
    support_audit: SupportLabelRenderAudit
    clades: CladeTableReport
    branch_lengths: BranchLengthDistributionReport
    support_rows: list[TreeSupportRow]
    branch_stats: TreeBranchStatisticsRow
    method_tier: MethodTierAssessment
    reviewer_summary: list[str]
    limitations: list[str]
    methods_summary: TreeValidationMethodsSummaryTextResult
    reviewer_audit_checklist: ReviewerAuditChecklist
    machine_manifest: dict[str, object]
