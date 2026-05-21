from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.discrete_mk import (
    DiscreteMkFitReport,
    compare_discrete_mk_model_ranking,
    fit_discrete_mk_model,
)
from bijux_phylogenetics.comparative.discrete_mode_recovery.reference_payloads import (
    GEIGER_FITDISCRETE_RECOVERY_REFERENCE_PAYLOADS,
)
from bijux_phylogenetics.simulation import (
    DiscreteHistoryRateRow,
    DiscreteTraitSimulationReport,
    simulate_discrete_histories,
)

from .comparison import _owner_model_name
from .formatting import _optional_string
from .models import (
    RECOVERY_TRAIT_NAME,
    DiscreteModeRecoveryExecutionRow,
    DiscreteModeRecoveryModelChoiceRow,
    DiscreteModeRecoveryScenario,
    DiscreteModeRecoveryWarningRow,
    _DiscreteRecoveryFitSnapshot,
)


def _simulate_case(
    tree_path: Path,
    scenario: DiscreteModeRecoveryScenario,
) -> DiscreteTraitSimulationReport:
    collection = simulate_discrete_histories(
        tree_path,
        states=scenario.states,
        rate_rows=scenario.rate_rows,
        root_state=scenario.root_state,
        transform=scenario.transform,
        transform_parameter_value=scenario.transform_parameter_value,
        replicates=1,
        seed=scenario.seed,
    )
    return collection.simulations[0]


def _build_bijux_fit_snapshot(
    tree_path: Path,
    traits_path: Path,
    scenario: DiscreteModeRecoveryScenario,
) -> _DiscreteRecoveryFitSnapshot:
    report = fit_discrete_mk_model(
        tree_path,
        traits_path,
        trait=RECOVERY_TRAIT_NAME,
        model=scenario.generating_model,
        transform=scenario.transform,
        lambda_bounds=scenario.lambda_bounds,
        kappa_bounds=scenario.kappa_bounds,
        delta_bounds=scenario.delta_bounds,
        early_burst_bounds=scenario.early_burst_bounds,
    )
    return _fit_snapshot_from_discrete_report(
        case_id=scenario.case_id,
        engine="bijux",
        report=report,
    )


def _fit_snapshot_from_discrete_report(
    *,
    case_id: str,
    engine: str,
    report: DiscreteMkFitReport,
) -> _DiscreteRecoveryFitSnapshot:
    warning_rows: list[DiscreteModeRecoveryWarningRow] = []
    if report.transform_fit is not None:
        for warning in report.transform_fit.warnings:
            warning_rows.append(
                DiscreteModeRecoveryWarningRow(
                    case_id=case_id,
                    recovery_engine=engine,
                    fitted_model=report.model,
                    kind=warning.kind,
                    message=warning.message,
                )
            )
    if report.overparameterized:
        warning_rows.append(
            DiscreteModeRecoveryWarningRow(
                case_id=case_id,
                recovery_engine=engine,
                fitted_model=report.model,
                kind="overparameterized",
                message=(
                    "The discrete Mk fit is overparameterized relative to the analyzed taxon count."
                ),
            )
        )
    if not report.optimizer_diagnostics.converged:
        warning_rows.append(
            DiscreteModeRecoveryWarningRow(
                case_id=case_id,
                recovery_engine=engine,
                fitted_model=report.model,
                kind="optimizer_not_converged",
                message=(
                    "The discrete Mk optimizer did not converge and should be interpreted cautiously."
                ),
            )
        )
    if (
        report.optimizer_diagnostics.hit_lower_parameter_bound
        or report.optimizer_diagnostics.hit_upper_parameter_bound
    ):
        warning_rows.append(
            DiscreteModeRecoveryWarningRow(
                case_id=case_id,
                recovery_engine=engine,
                fitted_model=report.model,
                kind="optimizer_boundary",
                message="One or more fitted discrete Mk rates hit an optimizer bound.",
            )
        )
    if (
        report.baseline_comparison is not None
        and report.baseline_comparison.preferred_model_by_aic == "equal-rates"
        and report.model != "equal-rates"
    ):
        warning_rows.append(
            DiscreteModeRecoveryWarningRow(
                case_id=case_id,
                recovery_engine=engine,
                fitted_model=report.model,
                kind="equal_rates_preferred",
                message=(
                    "The equal-rates baseline remains preferred by AIC over the fitted discrete Mk model."
                ),
            )
        )
    return _DiscreteRecoveryFitSnapshot(
        engine=engine,
        fitted_model=report.model,
        fit_status="ok",
        failure_reason=None,
        selected_model=None,
        parameter_name=(
            None
            if report.transform_fit is None
            else report.transform_fit.parameter_name
        ),
        parameter_value=(
            None
            if report.transform_fit is None
            else report.transform_fit.parameter_value
        ),
        parameter_count=report.parameter_count,
        log_likelihood=report.log_likelihood,
        aic=report.aic,
        aicc=report.aicc,
        overparameterized=report.overparameterized,
        optimizer_name=report.optimizer_diagnostics.optimizer_name,
        converged=report.optimizer_diagnostics.converged,
        hit_lower_parameter_bound=(
            report.optimizer_diagnostics.hit_lower_parameter_bound
            or (
                False
                if report.transform_fit is None
                else report.transform_fit.hit_lower_parameter_boundary
            )
        ),
        hit_upper_parameter_bound=(
            report.optimizer_diagnostics.hit_upper_parameter_bound
            or (
                False
                if report.transform_fit is None
                else report.transform_fit.hit_upper_parameter_boundary
            )
        ),
        rate_rows=[
            DiscreteHistoryRateRow(
                source_state=row.source_state,
                target_state=row.target_state,
                rate=row.rate,
            )
            for row in report.transition_rate_rows
        ],
        warning_rows=warning_rows,
    )


def _build_bijux_model_choice_rows(
    tree_path: Path,
    traits_path: Path,
    scenario: DiscreteModeRecoveryScenario,
) -> tuple[list[DiscreteModeRecoveryModelChoiceRow], str]:
    report = compare_discrete_mk_model_ranking(
        tree_path,
        traits_path,
        trait=RECOVERY_TRAIT_NAME,
        models=scenario.candidate_models,
        transform=scenario.transform,
        lambda_bounds=scenario.lambda_bounds,
        kappa_bounds=scenario.kappa_bounds,
        delta_bounds=scenario.delta_bounds,
        early_burst_bounds=scenario.early_burst_bounds,
    )
    return (
        [
            DiscreteModeRecoveryModelChoiceRow(
                case_id=scenario.case_id,
                generating_model=scenario.generating_model,
                transform=scenario.transform or "",
                recovery_engine="bijux",
                expected_selected_model=scenario.expected_selected_model,
                model=row.model,
                parameter_count=row.parameter_count,
                log_likelihood=row.log_likelihood,
                aic=row.aic,
                aicc=row.aicc,
                overparameterized=(
                    row.model == report.better_model
                    and "overparameterized" in " ".join(report.warnings)
                ),
                selected=row.selected,
            )
            for row in report.rows
        ],
        report.better_model,
    )


def _build_geiger_fit_snapshot(
    scenario: DiscreteModeRecoveryScenario,
    summary: dict[str, object],
    rows: list[dict[str, object]],
) -> _DiscreteRecoveryFitSnapshot:
    parameter_count = int(summary["parameter_count"])
    taxon_count = int(summary["taxon_count"])
    warning_rows: list[DiscreteModeRecoveryWarningRow] = []
    if parameter_count >= taxon_count:
        warning_rows.append(
            DiscreteModeRecoveryWarningRow(
                case_id=scenario.case_id,
                recovery_engine="geiger",
                fitted_model=scenario.generating_model,
                kind="overparameterized",
                message=(
                    "The stored geiger discrete fit is overparameterized relative to the analyzed taxon count."
                ),
            )
        )
    optimizer_result = summary.get("optimizer_result")
    optimizer_name = None
    converged = None
    if isinstance(optimizer_result, dict):
        optimizer_name = _optional_string(optimizer_result.get("best_method"))
        convergence_code = optimizer_result.get("convergence_code")
        if isinstance(convergence_code, int):
            converged = convergence_code == 0
        if converged is False:
            warning_rows.append(
                DiscreteModeRecoveryWarningRow(
                    case_id=scenario.case_id,
                    recovery_engine="geiger",
                    fitted_model=scenario.generating_model,
                    kind="optimizer_not_converged",
                    message="The stored geiger optimizer did not converge cleanly.",
                )
            )
    return _DiscreteRecoveryFitSnapshot(
        engine="geiger",
        fitted_model=scenario.generating_model,
        fit_status="ok",
        failure_reason=None,
        selected_model=None,
        parameter_name=_optional_string(summary.get("parameter_name")),
        parameter_value=(
            None
            if summary.get("parameter_value") is None
            else float(summary["parameter_value"])
        ),
        parameter_count=parameter_count,
        log_likelihood=float(summary["log_likelihood"]),
        aic=float(summary["aic"]),
        aicc=float(summary["aicc"]),
        overparameterized=parameter_count >= taxon_count,
        optimizer_name=optimizer_name,
        converged=converged,
        hit_lower_parameter_bound=(
            None
            if summary.get("hit_lower_parameter_boundary") is None
            else bool(summary["hit_lower_parameter_boundary"])
        ),
        hit_upper_parameter_bound=(
            None
            if summary.get("hit_upper_parameter_boundary") is None
            else bool(summary["hit_upper_parameter_boundary"])
        ),
        rate_rows=[
            DiscreteHistoryRateRow(
                source_state=str(row["source_state"]),
                target_state=str(row["target_state"]),
                rate=float(row["rate"]),
            )
            for row in rows
        ],
        warning_rows=warning_rows,
    )


def _build_geiger_model_choice_rows(
    scenario: DiscreteModeRecoveryScenario,
    summary: dict[str, object],
    rows: list[dict[str, object]] | dict[str, dict[str, object]],
) -> tuple[list[DiscreteModeRecoveryModelChoiceRow], str]:
    normalized_rows: list[dict[str, object]]
    if isinstance(rows, dict):
        normalized_rows = [dict(row) for row in rows.values()]
    else:
        normalized_rows = [dict(row) for row in rows]
    selection_rows = [
        DiscreteModeRecoveryModelChoiceRow(
            case_id=scenario.case_id,
            generating_model=scenario.generating_model,
            transform=scenario.transform or "",
            recovery_engine="geiger",
            expected_selected_model=scenario.expected_selected_model,
            model=_owner_model_name(str(row["model"])),
            parameter_count=int(row["parameter_count"]),
            log_likelihood=float(row["log_likelihood"]),
            aic=float(row["aic"]),
            aicc=float(row["aicc"]),
            overparameterized=int(row["parameter_count"])
            >= int(summary["taxon_count"]),
            selected=bool(row["selected"]),
        )
        for row in normalized_rows
    ]
    return selection_rows, _owner_model_name(str(summary["selected_model"]))


def _build_execution_row(
    case_id: str,
    snapshot: _DiscreteRecoveryFitSnapshot,
    *,
    operation: str,
) -> DiscreteModeRecoveryExecutionRow:
    return DiscreteModeRecoveryExecutionRow(
        case_id=case_id,
        recovery_engine=snapshot.engine,
        operation=operation,
        fitted_model=snapshot.fitted_model,
        fit_status=snapshot.fit_status,
        selected_model=snapshot.selected_model,
        optimizer_name=snapshot.optimizer_name,
        converged=snapshot.converged,
        hit_lower_parameter_bound=snapshot.hit_lower_parameter_bound,
        hit_upper_parameter_bound=snapshot.hit_upper_parameter_bound,
        overparameterized=snapshot.overparameterized,
        warning_count=len(snapshot.warning_rows),
        failure_reason=snapshot.failure_reason,
    )


def _build_model_comparison_execution_row(
    case_id: str,
    engine: str,
    selected_model: str,
    rows: list[DiscreteModeRecoveryModelChoiceRow],
) -> DiscreteModeRecoveryExecutionRow:
    return DiscreteModeRecoveryExecutionRow(
        case_id=case_id,
        recovery_engine=engine,
        operation="model-comparison",
        fitted_model="candidate-set",
        fit_status="ok",
        selected_model=selected_model,
        optimizer_name=None,
        converged=None,
        hit_lower_parameter_bound=None,
        hit_upper_parameter_bound=None,
        overparameterized=all(row.overparameterized for row in rows),
        warning_count=sum(1 for row in rows if row.overparameterized),
        failure_reason=None,
    )


def _geiger_payload_for_case(case_id: str) -> dict[str, object]:
    return GEIGER_FITDISCRETE_RECOVERY_REFERENCE_PAYLOADS[case_id]
