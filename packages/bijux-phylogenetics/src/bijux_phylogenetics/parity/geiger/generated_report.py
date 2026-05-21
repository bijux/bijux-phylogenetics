from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from bijux_phylogenetics.parity.geiger.generated_parity_report.contracts import (
    GeigerBenchmarkSummaryRow,
    GeigerBoundaryWarningSummaryRow,
    GeigerExcludedModelRow,
    GeigerGoalCoverageRow,
    GeigerOptimizerMismatchCategoryRow,
    GeigerSimulationRecoveryRow,
    GeigerToleranceRuleRow,
    GeneratedGeigerParityReport,
)
from bijux_phylogenetics.parity.geiger.generated_parity_report.governed_artifacts import (
    load_large_tree_benchmark_summary as _load_large_tree_benchmark_summary,
    load_real_dataset_benchmark_summary as _load_real_dataset_benchmark_summary,
    load_recovery_summary as _load_recovery_summary,
    load_sim_char_summary as _load_sim_char_summary,
    repository_root as _repository_root,
)
from bijux_phylogenetics.parity.geiger.generated_parity_report.coverage_policy import (
    boundary_warning_summaries as _boundary_warning_summaries,
    covered_models as _covered_models,
    excluded_models as _excluded_models,
    goal_coverage_rows as _goal_coverage_rows,
    optimizer_mismatch_categories as _optimizer_mismatch_categories,
    simulation_recovery_row as _simulation_recovery_row,
    tolerance_rules as _tolerance_rules,
)
from bijux_phylogenetics.parity.geiger.generated_parity_report.presentation import (
    write_generated_geiger_parity_report_json,
    write_generated_geiger_parity_report_markdown,
)

from .runner import GeigerParityReport, run_geiger_parity_cases

__all__ = [
    "GeneratedGeigerParityReport",
    "GeigerGoalCoverageRow",
    "GeigerExcludedModelRow",
    "GeigerToleranceRuleRow",
    "GeigerOptimizerMismatchCategoryRow",
    "GeigerBoundaryWarningSummaryRow",
    "GeigerSimulationRecoveryRow",
    "GeigerBenchmarkSummaryRow",
    "build_generated_geiger_parity_report",
    "write_generated_geiger_parity_report_json",
    "write_generated_geiger_parity_report_markdown",
]


def build_generated_geiger_parity_report(
    *,
    parity_report: GeigerParityReport | None = None,
    case_ids: list[str] | None = None,
    rscript_executable: str = "Rscript",
    failure_root: Path | None = None,
) -> GeneratedGeigerParityReport:
    """Build one generated geiger report from live parity and governed artifacts."""

    live_report = parity_report
    if live_report is None:
        live_report = run_geiger_parity_cases(
            case_ids=case_ids,
            rscript_executable=rscript_executable,
            failure_root=failure_root,
        )

    continuous_recovery = _load_recovery_summary(
        _repository_root()
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "resources"
        / "datasets"
        / "simulation"
        / "continuous_mode_recovery_panel"
        / "expected"
        / "workflow-summary.tsv"
    )
    discrete_recovery = _load_recovery_summary(
        _repository_root()
        / "packages"
        / "bijux-phylogenetics"
        / "src"
        / "bijux_phylogenetics"
        / "resources"
        / "datasets"
        / "simulation"
        / "discrete_mode_recovery_panel"
        / "expected"
        / "workflow-summary.tsv"
    )
    sim_char_report = _load_sim_char_summary()
    large_tree_summary = _load_large_tree_benchmark_summary()
    real_dataset_summary = _load_real_dataset_benchmark_summary()

    return GeneratedGeigerParityReport(
        generated_at_utc=datetime.now(UTC).isoformat(),
        goal_start=251,
        goal_end=289,
        r_version=_first_nonempty_string(
            observation.r_version for observation in live_report.observations
        ),
        geiger_version=_first_nonempty_string(
            observation.geiger_version for observation in live_report.observations
        ),
        live_case_count=live_report.case_count,
        live_passed_case_count=live_report.passed_case_count,
        live_failed_case_count=live_report.failed_case_count,
        live_skipped_case_count=live_report.skipped_case_count,
        all_live_cases_passed=live_report.all_passed,
        live_function_summary_rows=[
            {
                "function_name": row.function_name,
                "case_count": row.case_count,
                "passed_case_count": row.passed_case_count,
                "failed_case_count": row.failed_case_count,
                "skipped_case_count": row.skipped_case_count,
            }
            for row in live_report.summary_rows
        ],
        covered_models=_covered_models(),
        excluded_models=_excluded_models(),
        optimizer_mismatch_categories=_optimizer_mismatch_categories(live_report),
        tolerance_rules=_tolerance_rules(),
        boundary_warning_summaries=_boundary_warning_summaries(live_report),
        simulation_recovery_rows=[
            _simulation_recovery_row(
                continuous_recovery,
                panel_id="continuous_mode_recovery_panel",
                governed_value_success_key="parameter_pass_count",
                governed_value_row_key="parameter_row_count",
                governed_comparison_row_key="parameter_comparison_row_count",
                notes=[
                    "Stored governed expected bundle; generated from checked recovery panel outputs rather than rerunning the long simulation workflow."
                ],
            ),
            _simulation_recovery_row(
                discrete_recovery,
                panel_id="discrete_mode_recovery_panel",
                governed_value_success_key="rate_pass_count",
                governed_value_row_key="governed_rate_row_count",
                governed_comparison_row_key="governed_rate_comparison_row_count",
                notes=[
                    "Stored governed expected bundle; review-only ARD cases are retained instead of being misreported as successful recovery."
                ],
            ),
        ],
        benchmark_rows=[
            GeigerBenchmarkSummaryRow(
                benchmark_id="large_tree_model_fitting",
                case_count=large_tree_summary["case_count"],
                matched_case_count=large_tree_summary["geiger_match_case_count"],
                threshold_pass_case_count=large_tree_summary[
                    "threshold_pass_case_count"
                ],
                unstable_review_count=large_tree_summary["unstable_case_count"],
                too_slow_review_count=large_tree_summary["too_slow_case_count"],
                alignment_review_row_count=None,
                parity_row_count=None,
                notes=[
                    "Includes a real small-tier 100-taxon fit surface and preserves the heavy-tier timeout-review path rather than implying blanket large-tree success."
                ],
            ),
            GeigerBenchmarkSummaryRow(
                benchmark_id="real_dataset_macroevolution",
                case_count=real_dataset_summary["model_row_count"],
                matched_case_count=real_dataset_summary["selection_match_count"],
                threshold_pass_case_count=None,
                unstable_review_count=real_dataset_summary["unstable_review_count"],
                too_slow_review_count=0,
                alignment_review_row_count=real_dataset_summary[
                    "alignment_review_row_count"
                ],
                parity_row_count=real_dataset_summary["parity_row_count"],
                notes=[
                    "Real published dataset benchmark keeps sparse-state and near-boundary review lanes explicit instead of claiming universal row-level parity."
                ],
            ),
        ],
        sim_char_case_count=sim_char_report["case_count"],
        sim_char_all_passed=sim_char_report["all_passed"],
        goal_coverage_rows=_goal_coverage_rows(
            live_report=live_report,
            continuous_recovery=continuous_recovery,
            discrete_recovery=discrete_recovery,
            sim_char_report=sim_char_report,
            large_tree_summary=large_tree_summary,
            real_dataset_summary=real_dataset_summary,
        ),
        limitations=[
            *live_report.limitations,
            "This generated report counts live pass/fail/skip only for the consolidated geiger harness; dedicated stored-reference surfaces such as treedata, name.check, dtt, disparity, and rescale are inventoried separately and are not folded into the live totals.",
            "Simulation-recovery sections consume governed expected bundles so the report stays reproducible without rerunning the long recovery workflows on every invocation.",
        ],
    )
def _first_nonempty_string(values) -> str | None:
    for value in values:
        if isinstance(value, str) and value:
            return value
    return None
