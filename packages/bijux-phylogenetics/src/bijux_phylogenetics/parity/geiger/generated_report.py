from __future__ import annotations

from dataclasses import asdict, is_dataclass
from datetime import UTC, datetime
import json
from pathlib import Path

from bijux_phylogenetics.parity.geiger.generated_parity_report.contracts import (
    GeigerBenchmarkSummaryRow,
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

from .runner import GeigerParityReport, run_geiger_parity_cases


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


def write_generated_geiger_parity_report_json(
    path: Path,
    report: GeneratedGeigerParityReport,
) -> Path:
    """Write the generated geiger parity report as JSON."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_json_ready(report), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return path


def write_generated_geiger_parity_report_markdown(
    path: Path,
    report: GeneratedGeigerParityReport,
) -> Path:
    """Write the generated geiger parity report as Markdown."""

    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Geiger Parity Report",
        "",
        f"- Generated at: `{report.generated_at_utc}`",
        f"- Goal tranche: `{report.goal_start}`-`{report.goal_end}`",
        f"- R version: `{report.r_version or 'unavailable'}`",
        f"- geiger version: `{report.geiger_version or 'unavailable'}`",
        "",
        "## Live Summary",
        "",
        f"- Pass: `{report.live_passed_case_count}`",
        f"- Fail: `{report.live_failed_case_count}`",
        f"- Skip: `{report.live_skipped_case_count}`",
        f"- Total: `{report.live_case_count}`",
        f"- All passed: `{report.all_live_cases_passed}`",
        "",
        "### Live Function Coverage",
        "",
        "| Function | Cases | Passed | Failed | Skipped |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in report.live_function_summary_rows:
        lines.append(
            f"| `{row['function_name']}` | {row['case_count']} | {row['passed_case_count']} | {row['failed_case_count']} | {row['skipped_case_count']} |"
        )

    lines.extend(
        [
            "",
            "## Covered Models",
            "",
        ]
    )
    for model in report.covered_models:
        lines.append(f"- `{model}`")

    lines.extend(
        [
            "",
            "## Excluded Models",
            "",
            "| Goal | Surface | Exclusion Code | Reason |",
            "| --- | --- | --- | --- |",
        ]
    )
    for row in report.excluded_models:
        lines.append(
            f"| {row.goal_id} | `{row.surface}` | `{row.exclusion_code}` | {row.reason} |"
        )

    lines.extend(
        [
            "",
            "## Optimizer Mismatch Categories",
            "",
            "| Mismatch Type | Case Count | Case IDs |",
            "| --- | ---: | --- |",
        ]
    )
    for row in report.optimizer_mismatch_categories:
        lines.append(
            f"| `{row.mismatch_type}` | {row.case_count} | `{', '.join(row.case_ids)}` |"
        )

    lines.extend(
        [
            "",
            "## Tolerance Rules",
            "",
            "| Surface | Cases | Summary Tolerance | Row Policy | Field Overrides | Row Overrides |",
            "| --- | ---: | ---: | --- | --- | --- |",
        ]
    )
    for row in report.tolerance_rules:
        lines.append(
            "| "
            f"`{row.surface}` | {row.case_count} | {row.summary_tolerance} | "
            f"`{row.row_comparison_policy}` | `{json.dumps(row.field_tolerance_overrides, sort_keys=True)}` | "
            f"`{json.dumps(row.row_field_tolerance_overrides, sort_keys=True)}` |"
        )

    lines.extend(
        [
            "",
            "## Model-Boundary Warnings",
            "",
            "| Warning Kind | Case Count | Case IDs |",
            "| --- | ---: | --- |",
        ]
    )
    for row in report.boundary_warning_summaries:
        lines.append(
            f"| `{row.warning_kind}` | {row.case_count} | `{', '.join(row.case_ids)}` |"
        )

    lines.extend(
        [
            "",
            "## Simulation Recovery",
            "",
            "| Panel | Cases | Selection Review | Bijux Matches | geiger Matches | Governed Passes | Governed Rows | Comparison Rows | Expected Warnings | Present Warnings |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in report.simulation_recovery_rows:
        lines.append(
            f"| `{row.panel_id}` | {row.case_count} | {row.selection_review_case_count} | "
            f"{row.bijux_selection_match_count} | {row.geiger_selection_match_count} | "
            f"{row.governed_value_pass_count} | {row.governed_value_row_count} | "
            f"{row.governed_comparison_row_count} | {row.expected_warning_case_count} | "
            f"{row.expected_warning_present_count} |"
        )

    lines.extend(
        [
            "",
            "## Benchmarks",
            "",
            "| Benchmark | Cases | Matches | Threshold Passes | Unstable Reviews | Too Slow Reviews | Alignment Reviews | Parity Rows |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in report.benchmark_rows:
        lines.append(
            f"| `{row.benchmark_id}` | {row.case_count} | {row.matched_case_count} | "
            f"{'' if row.threshold_pass_case_count is None else row.threshold_pass_case_count} | "
            f"{row.unstable_review_count} | {row.too_slow_review_count} | "
            f"{'' if row.alignment_review_row_count is None else row.alignment_review_row_count} | "
            f"{'' if row.parity_row_count is None else row.parity_row_count} |"
        )

    lines.extend(
        [
            "",
            "## Goal Coverage",
            "",
            "| Goal | Surface | Status | Evidence | Cases | Passed | Failed | Skipped | Notes |",
            "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in report.goal_coverage_rows:
        lines.append(
            f"| {row.goal_id} | `{row.surface}` | `{row.status}` | `{row.evidence_kind}` | "
            f"{'' if row.observed_case_count is None else row.observed_case_count} | "
            f"{'' if row.passed_case_count is None else row.passed_case_count} | "
            f"{'' if row.failed_case_count is None else row.failed_case_count} | "
            f"{'' if row.skipped_case_count is None else row.skipped_case_count} | "
            f"{' '.join(row.notes)} |"
        )

    lines.extend(
        [
            "",
            "## Sim Char Envelope",
            "",
            f"- Cases: `{report.sim_char_case_count}`",
            f"- All passed: `{report.sim_char_all_passed}`",
            "",
            "## Limitations",
            "",
        ]
    )
    for item in report.limitations:
        lines.append(f"- {item}")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
def _first_nonempty_string(values) -> str | None:
    for value in values:
        if isinstance(value, str) and value:
            return value
    return None


def _json_ready(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        return {key: _json_ready(item) for key, item in asdict(value).items()}
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    return value
