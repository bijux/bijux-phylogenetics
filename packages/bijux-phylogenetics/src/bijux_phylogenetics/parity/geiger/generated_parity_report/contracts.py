from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class GeigerGoalCoverageRow:
    """One generated coverage row for the governed geiger goal tranche."""

    goal_id: int
    surface: str
    status: str
    evidence_kind: str
    observed_case_count: int | None
    passed_case_count: int | None
    failed_case_count: int | None
    skipped_case_count: int | None
    notes: list[str]


@dataclass(frozen=True, slots=True)
class GeigerExcludedModelRow:
    """One explicit non-coverage row in the governed geiger report."""

    goal_id: int
    surface: str
    exclusion_code: str
    reason: str


@dataclass(frozen=True, slots=True)
class GeigerToleranceRuleRow:
    """One distinct tolerance rule used by the live geiger parity registry."""

    surface: str
    case_count: int
    summary_tolerance: float
    row_comparison_policy: str
    field_tolerance_overrides: dict[str, float]
    row_field_tolerance_overrides: dict[str, float]
    notes: list[str]


@dataclass(frozen=True, slots=True)
class GeigerOptimizerMismatchCategoryRow:
    """One grouped optimizer mismatch category from the live geiger run."""

    mismatch_type: str
    case_count: int
    case_ids: list[str]


@dataclass(frozen=True, slots=True)
class GeigerBoundaryWarningSummaryRow:
    """One grouped boundary-warning kind summary from the live geiger run."""

    warning_kind: str
    case_count: int
    case_ids: list[str]


@dataclass(frozen=True, slots=True)
class GeigerSimulationRecoveryRow:
    """One generated simulation-recovery summary row."""

    panel_id: str
    case_count: int
    selection_review_case_count: int
    bijux_selection_match_count: int
    geiger_selection_match_count: int
    governed_value_pass_count: int
    governed_value_row_count: int
    governed_comparison_row_count: int
    expected_warning_case_count: int
    expected_warning_present_count: int
    notes: list[str]


@dataclass(frozen=True, slots=True)
class GeigerBenchmarkSummaryRow:
    """One generated benchmark summary row."""

    benchmark_id: str
    case_count: int
    matched_case_count: int
    threshold_pass_case_count: int | None
    unstable_review_count: int
    too_slow_review_count: int
    alignment_review_row_count: int | None
    parity_row_count: int | None
    notes: list[str]


@dataclass(slots=True)
class GeneratedGeigerParityReport:
    """Consolidated generated report over the governed geiger parity tranche."""

    generated_at_utc: str
    goal_start: int
    goal_end: int
    r_version: str | None
    geiger_version: str | None
    live_case_count: int
    live_passed_case_count: int
    live_failed_case_count: int
    live_skipped_case_count: int
    all_live_cases_passed: bool
    live_function_summary_rows: list[dict[str, object]]
    covered_models: list[str]
    excluded_models: list[GeigerExcludedModelRow]
    optimizer_mismatch_categories: list[GeigerOptimizerMismatchCategoryRow]
    tolerance_rules: list[GeigerToleranceRuleRow]
    boundary_warning_summaries: list[GeigerBoundaryWarningSummaryRow]
    simulation_recovery_rows: list[GeigerSimulationRecoveryRow]
    benchmark_rows: list[GeigerBenchmarkSummaryRow]
    sim_char_case_count: int
    sim_char_all_passed: bool
    goal_coverage_rows: list[GeigerGoalCoverageRow]
    limitations: list[str]
