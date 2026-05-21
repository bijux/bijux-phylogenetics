from __future__ import annotations

from dataclasses import asdict, is_dataclass
import json
from pathlib import Path

from .contracts import GeneratedGeigerParityReport


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
