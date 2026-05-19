from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime, timezone
import json
from pathlib import Path

from bijux_phylogenetics.benchmark import (
    benchmark_large_tree_model_fitting,
    benchmark_real_dataset_macroevolution,
)
from bijux_phylogenetics.parity.geiger.registry import list_geiger_parity_cases
from bijux_phylogenetics.simulation import validate_geiger_sim_char_reference_examples

from .runner import GeigerParityReport, run_geiger_parity_cases


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
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
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
                governed_value_pass_key="parameter_pass_count",
                governed_value_row_key="parameter_row_count",
                governed_comparison_row_key="parameter_comparison_row_count",
                notes=[
                    "Stored governed expected bundle; generated from checked recovery panel outputs rather than rerunning the long simulation workflow."
                ],
            ),
            _simulation_recovery_row(
                discrete_recovery,
                panel_id="discrete_mode_recovery_panel",
                governed_value_pass_key="rate_pass_count",
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


def _covered_models() -> list[str]:
    live_models = sorted(
        {case.function_name for case in list_geiger_parity_cases()},
    )
    extra_models = [
        "geiger::rescale(model='EB')",
        "geiger::rescale(model='delta')",
        "geiger::rescale(model='kappa')",
        "geiger::rescale(model='lambda')",
        "geiger::rescale(model='white')",
        "geiger::treedata",
        "geiger::name.check",
        "geiger::dtt",
        "geiger::disparity",
        "geiger::sim.char(model='BM')",
        "geiger::sim.char(model='discrete')",
        "geiger::sim.char(model='speciational')",
    ]
    return [*live_models, *extra_models]


def _excluded_models() -> list[GeigerExcludedModelRow]:
    return [
        GeigerExcludedModelRow(
            goal_id=261,
            surface="geiger::fitContinuous(trend aliases)",
            exclusion_code="fitcontinuous-trend-explicitly-excluded-this-round",
            reason=(
                "geiger exposes distinct rate_trend and mean_trend likelihoods, so Bijux does not claim one unambiguous trend parity contract."
            ),
        ),
        GeigerExcludedModelRow(
            goal_id=262,
            surface="geiger::fitContinuous standard errors",
            exclusion_code="fitcontinuous-standard-error-explicitly-excluded-this-round",
            reason=(
                "Bijux does not claim fitContinuous standard-error parity in this tranche and records that policy explicitly in the owned fit summaries."
            ),
        ),
        GeigerExcludedModelRow(
            goal_id=268,
            surface="geiger::fitDiscrete(model='meristic')",
            exclusion_code="fitdiscrete-meristic-explicitly-excluded-this-round",
            reason=(
                "Local geiger uses a distinct integer-state meristic contract and Bijux ordered-state Mk support is not claimed as meristic parity."
            ),
        ),
        GeigerExcludedModelRow(
            goal_id=278,
            surface="geiger::medusa",
            exclusion_code="geiger_medusa_explicitly_excluded_this_round",
            reason=(
                "Bijux does not yet implement MEDUSA stepwise rate-shift search, shift-count model growth, or branch-placement ranking."
            ),
        ),
        GeigerExcludedModelRow(
            goal_id=279,
            surface="geiger::bd.ms",
            exclusion_code="geiger_birth_death_explicitly_excluded_this_round",
            reason=(
                "The current owned diversification layer is a heuristic summary surface rather than a geiger::bd.ms-matched birth-death likelihood contract."
            ),
        ),
    ]


def _optimizer_mismatch_categories(
    report: GeigerParityReport,
) -> list[GeigerOptimizerMismatchCategoryRow]:
    rows_by_type: dict[str, list[str]] = {}
    for row in report.optimizer_triage_rows:
        rows_by_type.setdefault(row.mismatch_type, []).append(row.case_id)
    return [
        GeigerOptimizerMismatchCategoryRow(
            mismatch_type=mismatch_type,
            case_count=len(case_ids),
            case_ids=sorted(case_ids),
        )
        for mismatch_type, case_ids in sorted(rows_by_type.items())
    ]


def _tolerance_rules() -> list[GeigerToleranceRuleRow]:
    grouped: dict[
        tuple[str, float, str, str, str],
        list[str],
    ] = {}
    for case in list_geiger_parity_cases():
        key = (
            case.function_name,
            case.tolerance,
            case.row_comparison_policy,
            json.dumps(case.field_tolerances or {}, sort_keys=True),
            json.dumps(case.row_field_tolerances or {}, sort_keys=True),
        )
        grouped.setdefault(key, []).append(case.case_id)
    rows: list[GeigerToleranceRuleRow] = []
    for key, case_ids in sorted(grouped.items()):
        surface, tolerance, row_policy, field_overrides_json, row_overrides_json = key
        notes: list[str] = []
        if row_policy == "summary-only":
            notes.append("Rows are not claimed as row-level parity for this surface.")
        if field_overrides_json != "{}":
            notes.append("Summary fields use governed per-field tolerances.")
        if row_overrides_json != "{}":
            notes.append("Row fields use governed per-field tolerances.")
        rows.append(
            GeigerToleranceRuleRow(
                surface=surface,
                case_count=len(case_ids),
                summary_tolerance=tolerance,
                row_comparison_policy=row_policy,
                field_tolerance_overrides=json.loads(field_overrides_json),
                row_field_tolerance_overrides=json.loads(row_overrides_json),
                notes=notes,
            )
        )
    return rows


def _boundary_warning_summaries(
    report: GeigerParityReport,
) -> list[GeigerBoundaryWarningSummaryRow]:
    case_ids_by_kind: dict[str, set[str]] = {}
    for row in report.boundary_warning_rows:
        for kind in set(
            [*row.reference_boundary_warning_kinds, *row.bijux_boundary_warning_kinds]
        ):
            case_ids_by_kind.setdefault(kind, set()).add(row.case_id)
    return [
        GeigerBoundaryWarningSummaryRow(
            warning_kind=warning_kind,
            case_count=len(case_ids),
            case_ids=sorted(case_ids),
        )
        for warning_kind, case_ids in sorted(case_ids_by_kind.items())
    ]


def _simulation_recovery_row(
    payload: dict[str, int],
    *,
    panel_id: str,
    governed_value_pass_key: str,
    governed_value_row_key: str,
    governed_comparison_row_key: str,
    notes: list[str],
) -> GeigerSimulationRecoveryRow:
    return GeigerSimulationRecoveryRow(
        panel_id=panel_id,
        case_count=payload["case_count"],
        selection_review_case_count=payload["selection_review_case_count"],
        bijux_selection_match_count=payload["selection_match_count"],
        geiger_selection_match_count=payload["geiger_selection_match_count"],
        governed_value_pass_count=payload[governed_value_pass_key],
        governed_value_row_count=payload[governed_value_row_key],
        governed_comparison_row_count=payload[governed_comparison_row_key],
        expected_warning_case_count=payload["expected_warning_case_count"],
        expected_warning_present_count=payload["expected_warning_present_count"],
        notes=notes,
    )


def _goal_coverage_rows(
    *,
    live_report: GeigerParityReport,
    continuous_recovery: dict[str, int],
    discrete_recovery: dict[str, int],
    sim_char_report: dict[str, object],
    large_tree_summary: dict[str, int],
    real_dataset_summary: dict[str, int],
) -> list[GeigerGoalCoverageRow]:
    summary_by_function = {
        row.function_name: row for row in live_report.summary_rows
    }
    live_cases = list_geiger_parity_cases()
    continuous_fixture_count = len(
        {
            case.fixture_id
            for case in live_cases
            if case.function_name.startswith("geiger::fitContinuous")
        }
    )
    discrete_fixture_count = len(
        {
            case.fixture_id
            for case in live_cases
            if case.function_name.startswith("geiger::fitDiscrete")
        }
    )

    def live_row(function_name: str, goal_id: int) -> GeigerGoalCoverageRow:
        row = summary_by_function[function_name]
        return GeigerGoalCoverageRow(
            goal_id=goal_id,
            surface=function_name,
            status="live-parity",
            evidence_kind="live-geiger-harness",
            observed_case_count=row.case_count,
            passed_case_count=row.passed_case_count,
            failed_case_count=row.failed_case_count,
            skipped_case_count=row.skipped_case_count,
            notes=[],
        )

    def aggregate_live_row(
        *,
        goal_id: int,
        surface: str,
        function_names: list[str],
        notes: list[str] | None = None,
    ) -> GeigerGoalCoverageRow:
        matched_rows = [summary_by_function[name] for name in function_names]
        return GeigerGoalCoverageRow(
            goal_id=goal_id,
            surface=surface,
            status="live-parity",
            evidence_kind="live-geiger-harness",
            observed_case_count=sum(row.case_count for row in matched_rows),
            passed_case_count=sum(row.passed_case_count for row in matched_rows),
            failed_case_count=sum(row.failed_case_count for row in matched_rows),
            skipped_case_count=sum(row.skipped_case_count for row in matched_rows),
            notes=[] if notes is None else notes,
        )

    supporting_transform_case_count = sum(
        row.case_count
        for name, row in summary_by_function.items()
        if "model='lambda'" in name
        or "model='kappa'" in name
        or "model='delta'" in name
        or "model='EB'" in name
        or "model='white'" in name
        or "transform='lambda'" in name
        or "transform='kappa'" in name
        or "transform='delta'" in name
        or "transform='EB'" in name
    )

    return [
        GeigerGoalCoverageRow(
            goal_id=251,
            surface="live geiger parity harness",
            status="live-parity",
            evidence_kind="live-geiger-harness",
            observed_case_count=live_report.case_count,
            passed_case_count=live_report.passed_case_count,
            failed_case_count=live_report.failed_case_count,
            skipped_case_count=live_report.skipped_case_count,
            notes=[],
        ),
        GeigerGoalCoverageRow(
            goal_id=252,
            surface="shared continuous fixture catalog",
            status="artifact-backed",
            evidence_kind="live-fixture-inventory",
            observed_case_count=continuous_fixture_count,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Count reflects distinct governed continuous fixture ids used by the live geiger harness."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=253,
            surface="shared discrete fixture catalog",
            status="artifact-backed",
            evidence_kind="live-fixture-inventory",
            observed_case_count=discrete_fixture_count,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Count reflects distinct governed discrete fixture ids used by the live geiger harness."
            ],
        ),
        live_row("geiger::fitContinuous(model='BM')", 254),
        live_row("geiger::fitContinuous(model='OU')", 255),
        live_row("geiger::fitContinuous(model='EB')", 256),
        live_row("geiger::fitContinuous(model='lambda')", 257),
        live_row("geiger::fitContinuous(model='kappa')", 258),
        live_row("geiger::fitContinuous(model='delta')", 259),
        live_row("geiger::fitContinuous(model='white')", 260),
        GeigerGoalCoverageRow(
            goal_id=261,
            surface="geiger::fitContinuous(trend aliases)",
            status="explicit-exclusion",
            evidence_kind="governed-exclusion",
            observed_case_count=None,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Distinct rate_trend and mean_trend likelihoods are not collapsed into one claimed parity lane."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=262,
            surface="geiger::fitContinuous standard errors",
            status="explicit-exclusion",
            evidence_kind="governed-exclusion",
            observed_case_count=None,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Standard-error parity remains excluded and is surfaced explicitly in owned continuous fit outputs."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=263,
            surface="fitContinuous bounds and control review",
            status="live-parity",
            evidence_kind="live-geiger-harness",
            observed_case_count=8,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Governed lower-bound, upper-bound, weak-signal, and control-review cases are included in the live fitContinuous suite."
            ],
        ),
        live_row("geiger::fitContinuous(model comparison)", 264),
        live_row("geiger::fitDiscrete(model='ER')", 265),
        live_row("geiger::fitDiscrete(model='SYM')", 266),
        live_row("geiger::fitDiscrete(model='ARD')", 267),
        GeigerGoalCoverageRow(
            goal_id=268,
            surface="geiger::fitDiscrete(model='meristic')",
            status="explicit-exclusion",
            evidence_kind="governed-exclusion",
            observed_case_count=None,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Meristic parity is excluded rather than being mislabeled as generic ordered-state Mk."
            ],
        ),
        live_row("geiger::fitDiscrete(model='ER', transform='lambda')", 269),
        aggregate_live_row(
            goal_id=270,
            surface="geiger::fitDiscrete(transform='kappa')",
            function_names=[
                "geiger::fitDiscrete(model='ER', transform='kappa')",
                "geiger::fitDiscrete(model='SYM', transform='kappa')",
            ],
        ),
        aggregate_live_row(
            goal_id=271,
            surface="geiger::fitDiscrete(transform='delta')",
            function_names=[
                "geiger::fitDiscrete(model='ER', transform='delta')",
                "geiger::fitDiscrete(model='SYM', transform='delta')",
            ],
        ),
        live_row("geiger::fitDiscrete(model='ER', transform='EB')", 272),
        GeigerGoalCoverageRow(
            goal_id=273,
            surface="geiger::rescale",
            status="artifact-backed",
            evidence_kind="shared-transform-support",
            observed_case_count=supporting_transform_case_count,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "The shared rescale engine is exercised by the governed transform fit lanes that feed this live suite and the dedicated rescale parity surface."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=274,
            surface="geiger::treedata",
            status="artifact-backed",
            evidence_kind="supporting-owner-surface",
            observed_case_count=live_report.case_count,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "The shared tree-trait alignment owner path is a prerequisite for all live fit cases, but this report does not fold its dedicated stored-reference parity rows into live pass/fail totals."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=275,
            surface="geiger::name.check",
            status="artifact-backed",
            evidence_kind="supporting-owner-surface",
            observed_case_count=real_dataset_summary["alignment_review_row_count"],
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Name-check evidence is surfaced through the shared alignment review path and dedicated stored-reference parity tests."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=276,
            surface="geiger::dtt",
            status="artifact-backed",
            evidence_kind="stored-reference-governed",
            observed_case_count=None,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Observed-curve parity is governed outside the consolidated live runner and is inventoried here as a dedicated stored-reference surface."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=277,
            surface="geiger::disparity",
            status="artifact-backed",
            evidence_kind="stored-reference-governed",
            observed_case_count=None,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Direct clade-disparity parity is governed outside the consolidated live runner and is inventoried here as a dedicated stored-reference surface."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=278,
            surface="geiger::medusa",
            status="explicit-exclusion",
            evidence_kind="governed-exclusion",
            observed_case_count=None,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "MEDUSA-style stepwise rate-shift search remains explicitly excluded."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=279,
            surface="geiger::bd.ms",
            status="explicit-exclusion",
            evidence_kind="governed-exclusion",
            observed_case_count=None,
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "The owned diversification summary layer is not claimed as bd.ms likelihood parity."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=280,
            surface="geiger::sim.char",
            status="artifact-backed",
            evidence_kind="simulation-envelope-validation",
            observed_case_count=int(sim_char_report["case_count"]),
            passed_case_count=int(sim_char_report["case_count"])
            if bool(sim_char_report["all_passed"])
            else None,
            failed_case_count=0 if bool(sim_char_report["all_passed"]) else None,
            skipped_case_count=0,
            notes=[
                "Envelope validation is summary-level rather than exact draw-level parity."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=281,
            surface="fitContinuous recovery benchmark",
            status="artifact-backed",
            evidence_kind="simulation-recovery-bundle",
            observed_case_count=continuous_recovery["case_count"],
            passed_case_count=continuous_recovery["selection_match_count"],
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Stored governed workflow bundle captures strong, weak-identification, and transformed review cases."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=282,
            surface="fitDiscrete recovery benchmark",
            status="artifact-backed",
            evidence_kind="simulation-recovery-bundle",
            observed_case_count=discrete_recovery["case_count"],
            passed_case_count=discrete_recovery["selection_match_count"],
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Stable ER and SYM lanes are separated from ARD weak-identification review surfaces."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=283,
            surface="optimizer disagreement triage",
            status="live-parity",
            evidence_kind="live-row-registry",
            observed_case_count=len(live_report.optimizer_triage_rows),
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[],
        ),
        GeigerGoalCoverageRow(
            goal_id=284,
            surface="parameterization-difference registry",
            status="live-parity",
            evidence_kind="live-row-registry",
            observed_case_count=len(live_report.parameterization_registry_rows),
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[],
        ),
        GeigerGoalCoverageRow(
            goal_id=285,
            surface="likelihood-constant policy",
            status="live-parity",
            evidence_kind="live-row-registry",
            observed_case_count=len(live_report.likelihood_policy_rows),
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[],
        ),
        GeigerGoalCoverageRow(
            goal_id=286,
            surface="model-boundary warnings",
            status="live-parity",
            evidence_kind="live-row-registry",
            observed_case_count=len(live_report.boundary_warning_rows),
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[],
        ),
        GeigerGoalCoverageRow(
            goal_id=287,
            surface="AIC weights and model confidence",
            status="live-parity",
            evidence_kind="live-row-registry",
            observed_case_count=len(live_report.model_confidence_rows),
            passed_case_count=None,
            failed_case_count=None,
            skipped_case_count=None,
            notes=[],
        ),
        GeigerGoalCoverageRow(
            goal_id=288,
            surface="large-tree model fitting benchmark",
            status="artifact-backed",
            evidence_kind="benchmark-surface",
            observed_case_count=large_tree_summary["case_count"],
            passed_case_count=large_tree_summary["geiger_match_case_count"],
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Heavy-tier timeout-review is retained as benchmark evidence instead of being hidden as a skipped case."
            ],
        ),
        GeigerGoalCoverageRow(
            goal_id=289,
            surface="real-dataset macroevolution benchmark",
            status="artifact-backed",
            evidence_kind="benchmark-surface",
            observed_case_count=real_dataset_summary["model_row_count"],
            passed_case_count=real_dataset_summary["selection_match_count"],
            failed_case_count=None,
            skipped_case_count=None,
            notes=[
                "Selection-level agreement is reported separately from sparse-state and near-boundary review rows."
            ],
        ),
    ]


def _load_recovery_summary(path: Path) -> dict[str, int | str]:
    rows = _read_tsv_rows(path)
    if len(rows) != 1:
        raise ValueError(f"expected one workflow summary row in {path}")
    row = rows[0]
    payload: dict[str, int | str] = {}
    for key, value in row.items():
        if key == "dataset_id":
            payload[key] = value
            continue
        payload[key] = int(value)
    return payload


def _load_sim_char_summary() -> dict[str, object]:
    report = validate_geiger_sim_char_reference_examples()
    return {
        "case_count": report.case_count,
        "all_passed": report.all_passed,
    }


def _load_large_tree_benchmark_summary() -> dict[str, int]:
    report = benchmark_large_tree_model_fitting(tier="small")
    return {
        "case_count": report.case_count,
        "geiger_match_case_count": report.geiger_match_case_count,
        "threshold_pass_case_count": report.threshold_pass_case_count,
        "too_slow_case_count": report.too_slow_case_count,
        "unstable_case_count": report.unstable_case_count,
    }


def _load_real_dataset_benchmark_summary() -> dict[str, int]:
    report = benchmark_real_dataset_macroevolution()
    return {
        "summary_row_count": len(report.summary_rows),
        "model_row_count": len(report.model_rows),
        "alignment_review_row_count": len(report.alignment_review_rows),
        "parity_row_count": len(report.parity_rows),
        "selection_match_count": sum(
            1 for row in report.summary_rows if row.selection_matches_geiger
        ),
        "unstable_review_count": sum(
            1 for row in report.summary_rows if not row.stable_conclusion_supported
        ),
    }


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _first_nonempty_string(values) -> str | None:
    for value in values:
        if isinstance(value, str) and value:
            return value
    return None


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[6]


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
