from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.reporting import (
    ComparativeMethodReport,
    ComparativeMethodsSummaryTextResult,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    MethodTierAssessment,
)
from bijux_phylogenetics.reports.review import ReviewerAuditChecklist


@dataclass(frozen=True, slots=True)
class ComparativeAnalysisSummaryRow:
    response: str
    formula: str
    predictor_count: int
    analysis_taxa: int
    excluded_taxa: int
    selected_model: str
    pgls_lambda: float
    pgls_log_likelihood: float
    pgls_r_squared: float
    phylogenetic_signal_k: float
    phylogenetic_signal_lambda: float
    independent_contrast_count: int
    better_model_aicc_delta: float


@dataclass(frozen=True, slots=True)
class ComparativeCoefficientTableRow:
    term: str
    estimate: float
    standard_error: float
    test_statistic: float
    p_value: float
    lower_95_confidence_interval: float
    upper_95_confidence_interval: float
    degrees_of_freedom: int
    inference_distribution: str
    significant: bool


@dataclass(frozen=True, slots=True)
class ComparativeResidualTableRow:
    analysis: str
    residual_variance: float
    max_abs_standardized_residual: float
    phylogenetic_residual_lambda: float | None
    max_leverage: float | None
    outlier_taxa: tuple[str, ...]
    high_leverage_taxa: tuple[str, ...]
    warnings: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ComparativeSignalTableRow:
    trait: str
    taxon_count: int
    blombergs_k: float
    pagels_lambda: float
    lambda_log_likelihood: float
    lambda_null_log_likelihood: float
    lambda_brownian_log_likelihood: float
    independent_contrast_count: int
    independent_contrast_root_estimate: float


@dataclass(frozen=True, slots=True)
class ComparativeInterpretationRow:
    topic: str
    claim: str
    evidence: str
    caution: str


@dataclass(frozen=True, slots=True)
class ComparativeAuditTableRow:
    analysis: str
    taxa_used: tuple[str, ...]
    traits_used: tuple[str, ...]
    excluded_taxa: tuple[str, ...]
    assumptions: tuple[str, ...]
    warnings: tuple[str, ...]


@dataclass(slots=True)
class ComparativeReportPackageResult:
    output_dir: Path
    report_path: Path
    methods_summary_path: Path
    reviewer_audit_checklist_path: Path
    summary_table_path: Path
    coefficient_table_path: Path
    residual_table_path: Path
    signal_table_path: Path
    model_comparison_table_path: Path
    interpretation_table_path: Path
    audit_table_path: Path
    contrast_table_path: Path
    manifest_path: Path
    report: ComparativeMethodReport
    methods_summary: ComparativeMethodsSummaryTextResult
    summary_row: ComparativeAnalysisSummaryRow
    coefficient_rows: list[ComparativeCoefficientTableRow]
    residual_rows: list[ComparativeResidualTableRow]
    signal_row: ComparativeSignalTableRow
    interpretation_rows: list[ComparativeInterpretationRow]
    audit_rows: list[ComparativeAuditTableRow]
    method_tier: MethodTierAssessment
    reviewer_audit_checklist: ReviewerAuditChecklist
    machine_manifest: dict[str, object]
