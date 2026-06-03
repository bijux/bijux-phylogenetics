from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.bayesian.beast.models import (
    FossilCalibrationValidationReport,
    TipDatingValidationReport,
)
from bijux_phylogenetics.datasets.study_inputs import MetadataColumnCompleteness
from bijux_phylogenetics.phylo.alignment import AlignmentForensicReport


@dataclass(slots=True)
class DatasetReadinessSummary:
    """Readiness summary for a tree plus linked metadata and trait tables."""

    tree_path: Path
    metadata_path: Path
    traits_path: Path
    tree_taxa: int
    analysis_taxa: list[str]
    missing_metadata_taxa: list[str]
    missing_trait_taxa: list[str]
    metadata_only_taxa: list[str]
    trait_only_taxa: list[str]
    metadata_column_completeness: list[MetadataColumnCompleteness]
    unusable_trait_columns: list[str]
    ready_for_comparative_analysis: bool
    blockers: list[str]
    warnings: list[str]


@dataclass(slots=True)
class DatasetAuditFinding:
    """One categorized dataset blocker or warning."""

    severity: str
    category: str
    code: str
    message: str
    affected_analyses: list[str]


@dataclass(slots=True)
class DatasetAnalysisDecision:
    """Eligibility decision for one downstream analysis family."""

    analysis: str
    decision: str
    reasons: list[str]


@dataclass(slots=True)
class DatasetReadinessLevel:
    """Named dataset readiness level suitable for reviewer-facing summaries."""

    level: str
    decision: str
    reasons: list[str]


@dataclass(slots=True)
class DatasetCrosswalkRow:
    """Explicit mapping for one taxon across dataset surfaces."""

    taxon: str
    tree_tip: str | None
    alignment_id: str | None
    metadata_id: str | None
    trait_id: str | None
    tip_date_id: str | None
    geography_source: str | None
    calibration_targets: list[str]
    external_taxonomy_ids: dict[str, str]


@dataclass(slots=True)
class DatasetCrosswalkReport:
    """Explicit taxon crosswalk across tree, alignment, metadata, traits, and dates."""

    tree_path: Path
    metadata_path: Path
    traits_path: Path
    alignment_path: Path | None
    tip_dates_path: Path | None
    calibration_path: Path | None
    rows: list[DatasetCrosswalkRow]


@dataclass(slots=True)
class DatasetCompletenessRow:
    """Presence or absence of one taxon across dataset surfaces."""

    taxon: str
    in_tree: bool
    in_alignment: bool
    in_metadata: bool
    in_traits: bool
    in_tip_dates: bool
    in_geography: bool
    in_calibrations: bool


@dataclass(slots=True)
class DatasetCompletenessMatrix:
    """Taxon-by-surface completeness summary."""

    tree_path: Path
    metadata_path: Path
    traits_path: Path
    alignment_path: Path | None
    tip_dates_path: Path | None
    calibration_path: Path | None
    geography_columns: list[str]
    rows: list[DatasetCompletenessRow]
    surface_counts: dict[str, int]


@dataclass(slots=True)
class DatasetExclusionRow:
    """One excluded taxon with explicit causes and downstream consequences."""

    taxon: str
    causes: list[str]
    first_failed_surface: str
    affected_analyses: list[str]


@dataclass(slots=True)
class DatasetExclusionTable:
    """All excluded taxa with exact causes and affected analysis families."""

    rows: list[DatasetExclusionRow]


@dataclass(slots=True)
class DatasetOrderingConflict:
    """One ordering conflict against the canonical taxon order."""

    surface: str
    taxon: str
    expected_index: int
    observed_index: int


@dataclass(slots=True)
class DatasetOrderingAudit:
    """Ordering audit for tree-linked surfaces."""

    canonical_surface: str
    consistent: bool
    drifted_surfaces: list[str]
    conflicts: list[DatasetOrderingConflict]


@dataclass(slots=True)
class DatasetPruningStepSummary:
    """Sample-size summary after one dataset linkage or pruning step."""

    step: str
    input_taxa: int
    retained_taxa: int
    excluded_taxa: int
    reason: str


@dataclass(slots=True)
class DatasetGroupImbalanceWarning:
    """Warning when pruning removes most taxa from one categorical group."""

    surface: str
    group_column: str
    group: str
    original_count: int
    retained_count: int
    removed_count: int
    removed_fraction: float
    message: str


@dataclass(frozen=True, slots=True)
class DatasetMismatchRow:
    """One taxon with inconsistent presence across requested dataset surfaces."""

    taxon: str
    present_surfaces: list[str]
    missing_surfaces: list[str]
    message: str


@dataclass(slots=True)
class DatasetMismatchReport:
    """Taxon-by-surface mismatch report across the requested inputs."""

    requested_surfaces: list[str]
    rows: list[DatasetMismatchRow]
    mismatch_counts: dict[str, int]


@dataclass(frozen=True, slots=True)
class DatasetRiskComponent:
    """Transparent contribution of one surface to overall dataset risk."""

    component: str
    score: float
    reasons: list[str]


@dataclass(slots=True)
class DatasetRiskScoreReport:
    """Dataset risk score decomposed across major evidence surfaces."""

    total_score: float
    risk_level: str
    components: list[DatasetRiskComponent]


@dataclass(frozen=True, slots=True)
class DatasetFixRecommendation:
    """One minimal change that unlocks additional downstream analyses."""

    priority: int
    summary: str
    affected_surfaces: list[str]
    unlocks_analyses: list[str]


@dataclass(slots=True)
class DatasetMinimalFixPlan:
    """Smallest next fixes that materially improve dataset readiness."""

    recommendations: list[DatasetFixRecommendation]


@dataclass(frozen=True, slots=True)
class DatasetReviewerChecklistItem:
    """Reviewer-facing pass/risk/block item for one audit surface."""

    section: str
    status: str
    summary: str
    evidence: list[str]


@dataclass(slots=True)
class DatasetReviewerChecklist:
    """Reviewer-readable checklist across all major dataset trust surfaces."""

    items: list[DatasetReviewerChecklistItem]


@dataclass(slots=True)
class DatasetAuditReport:
    """Integrated audit across tree, alignment, metadata, traits, and optional dating surfaces."""

    tree_path: Path
    metadata_path: Path
    traits_path: Path
    alignment_path: Path | None
    tip_dates_path: Path | None
    calibration_path: Path | None
    tree_taxa: int
    analysis_taxa: list[str]
    readiness_decision: str
    allowed_analyses: list[str]
    blocked_analyses: list[str]
    blockers: list[str]
    warnings: list[str]
    blocker_categories: dict[str, list[str]]
    warning_categories: dict[str, list[str]]
    findings: list[DatasetAuditFinding]
    analysis_decisions: list[DatasetAnalysisDecision]
    readiness_levels: list[DatasetReadinessLevel]
    crosswalk: DatasetCrosswalkReport
    completeness_matrix: DatasetCompletenessMatrix
    exclusion_table: DatasetExclusionTable
    ordering_audit: DatasetOrderingAudit
    pruning_steps: list[DatasetPruningStepSummary]
    group_imbalance_warnings: list[DatasetGroupImbalanceWarning]
    dataset_readiness: DatasetReadinessSummary
    alignment_forensic: AlignmentForensicReport | None
    tip_dates: TipDatingValidationReport | None
    calibrations: FossilCalibrationValidationReport | None
    mismatch_report: DatasetMismatchReport
    risk_score: DatasetRiskScoreReport
    minimal_fix_plan: DatasetMinimalFixPlan
    reviewer_checklist: DatasetReviewerChecklist
