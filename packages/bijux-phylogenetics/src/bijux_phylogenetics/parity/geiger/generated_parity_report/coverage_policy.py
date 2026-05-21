from __future__ import annotations

import json

from bijux_phylogenetics.parity.geiger.registry import list_geiger_parity_cases
from bijux_phylogenetics.parity.geiger.runner import GeigerParityReport

from .contracts import (
    GeigerBoundaryWarningSummaryRow,
    GeigerExcludedModelRow,
    GeigerGoalCoverageRow,
    GeigerOptimizerMismatchCategoryRow,
    GeigerSimulationRecoveryRow,
    GeigerToleranceRuleRow,
)


def covered_models() -> list[str]:
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


def excluded_models() -> list[GeigerExcludedModelRow]:
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


def optimizer_mismatch_categories(
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


def tolerance_rules() -> list[GeigerToleranceRuleRow]:
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
        if "fitdiscrete-lambda-weak-signal-review" in case_ids:
            notes.append(
                "The weak lambda review compares the stable likelihood, AIC, AICc, and pruning surface rather than the nondeterministic plateau lambda value from live geiger."
            )
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


def boundary_warning_summaries(
    report: GeigerParityReport,
) -> list[GeigerBoundaryWarningSummaryRow]:
    case_ids_by_kind: dict[str, set[str]] = {}
    for row in report.boundary_warning_rows:
        for kind in {
            *row.reference_boundary_warning_kinds,
            *row.bijux_boundary_warning_kinds,
        }:
            case_ids_by_kind.setdefault(kind, set()).add(row.case_id)
    return [
        GeigerBoundaryWarningSummaryRow(
            warning_kind=warning_kind,
            case_count=len(case_ids),
            case_ids=sorted(case_ids),
        )
        for warning_kind, case_ids in sorted(case_ids_by_kind.items())
    ]


def simulation_recovery_row(
    payload: dict[str, int],
    *,
    panel_id: str,
    governed_value_success_key: str,
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
        governed_value_pass_count=payload[governed_value_success_key],
        governed_value_row_count=payload[governed_value_row_key],
        governed_comparison_row_count=payload[governed_comparison_row_key],
        expected_warning_case_count=payload["expected_warning_case_count"],
        expected_warning_present_count=payload["expected_warning_present_count"],
        notes=notes,
    )


def goal_coverage_rows(
    *,
    live_report: GeigerParityReport,
    continuous_recovery: dict[str, int],
    discrete_recovery: dict[str, int],
    sim_char_report: dict[str, object],
    large_tree_summary: dict[str, int],
    real_dataset_summary: dict[str, int],
) -> list[GeigerGoalCoverageRow]:
    summary_by_function = {row.function_name: row for row in live_report.summary_rows}
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

    def live_row(
        function_name: str,
        goal_id: int,
        *,
        notes: list[str] | None = None,
    ) -> GeigerGoalCoverageRow:
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
            notes=[] if notes is None else notes,
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
        live_row(
            "geiger::fitDiscrete(model='ER', transform='lambda')",
            269,
            notes=[
                "The weak-signal lambda review is governed by objective-surface agreement because live geiger returns nondeterministic plateau lambda values while Bijux resolves the same flat surface to the lower boundary."
            ],
        ),
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
