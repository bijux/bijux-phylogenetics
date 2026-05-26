from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.evolutionary_modes import (
    compare_fitcontinuous_model_ranking,
)

from .formatting import _optional_string
from .models import (
    RECOVERY_TRAIT_NAME,
    ContinuousModeRecoveryModelChoiceRow,
    ContinuousModeRecoveryParameterComparisonRow,
    ContinuousModeRecoveryParameterRow,
    ContinuousModeRecoveryScenario,
    _RecoveryComparisonSnapshot,
    _RecoveryFitSnapshot,
)


def _build_bijux_comparison_snapshot(
    tree_path: Path,
    traits_path: Path,
    scenario: ContinuousModeRecoveryScenario,
) -> _RecoveryComparisonSnapshot:
    report = compare_fitcontinuous_model_ranking(
        tree_path,
        traits_path,
        trait=RECOVERY_TRAIT_NAME,
        modes=scenario.candidate_modes,
        lambda_bounds=scenario.lambda_bounds,
        kappa_bounds=scenario.kappa_bounds,
        delta_bounds=scenario.delta_bounds,
        ou_bounds=scenario.ou_bounds,
        early_burst_bounds=scenario.early_burst_bounds,
    )
    rows = [
        ContinuousModeRecoveryModelChoiceRow(
            case_id=scenario.case_id,
            generating_model=scenario.generating_model,
            recovery_engine="bijux",
            expected_selected_model=scenario.expected_selected_model,
            model=row.model,
            parameter_count=row.parameter_count,
            log_likelihood=row.log_likelihood,
            aic=row.aic,
            aicc=row.aicc,
            selected=row.selected,
        )
        for row in report.rows
    ]
    return _RecoveryComparisonSnapshot(
        engine="bijux",
        fit_status="ok",
        failure_reason=None,
        selected_model=report.better_model,
        rows=rows,
        warning_count=len(report.warnings),
    )


def _build_geiger_comparison_snapshot(
    scenario: ContinuousModeRecoveryScenario,
    summary: dict[str, object],
    rows: list[dict[str, object]],
) -> _RecoveryComparisonSnapshot:
    return _RecoveryComparisonSnapshot(
        engine="geiger",
        fit_status="ok",
        failure_reason=None,
        selected_model=_optional_string(summary.get("selected_model")),
        rows=[
            ContinuousModeRecoveryModelChoiceRow(
                case_id=scenario.case_id,
                generating_model=scenario.generating_model,
                recovery_engine="geiger",
                expected_selected_model=scenario.expected_selected_model,
                model=str(row["model"]),
                parameter_count=int(row["parameter_count"]),
                log_likelihood=float(row["log_likelihood"]),
                aic=float(row["aic"]),
                aicc=float(row["aicc"]),
                selected=bool(row["selected"]),
            )
            for row in rows
        ],
        warning_count=int(summary.get("warning_count", 0)),
    )


def _build_parameter_rows(
    *,
    scenario: ContinuousModeRecoveryScenario,
    bijux_fit: _RecoveryFitSnapshot,
    geiger_fit: _RecoveryFitSnapshot,
) -> list[ContinuousModeRecoveryParameterRow]:
    rows: list[ContinuousModeRecoveryParameterRow] = []
    truth_by_parameter = _truth_parameters(scenario)
    for parameter, true_value in truth_by_parameter.items():
        tolerance = scenario.parameter_tolerances.get(parameter)
        if tolerance is None:
            continue
        bijux_estimate = _estimate_for_parameter(bijux_fit, parameter)
        if bijux_estimate is not None:
            rows.append(
                _parameter_row(
                    case_id=scenario.case_id,
                    generating_model=scenario.generating_model,
                    recovery_engine="bijux",
                    fitted_model=scenario.generating_model,
                    fit_status=bijux_fit.fit_status,
                    parameter=parameter,
                    true_value=true_value,
                    estimated_value=bijux_estimate,
                    tolerance=tolerance,
                )
            )
        geiger_estimate = _estimate_for_parameter(geiger_fit, parameter)
        if geiger_estimate is not None:
            rows.append(
                _parameter_row(
                    case_id=scenario.case_id,
                    generating_model=scenario.generating_model,
                    recovery_engine="geiger",
                    fitted_model=scenario.generating_model,
                    fit_status=geiger_fit.fit_status,
                    parameter=parameter,
                    true_value=true_value,
                    estimated_value=geiger_estimate,
                    tolerance=tolerance,
                )
            )
    return rows


def _build_parameter_comparison_rows(
    parameter_rows: list[ContinuousModeRecoveryParameterRow],
) -> list[ContinuousModeRecoveryParameterComparisonRow]:
    grouped: dict[tuple[str, str], list[ContinuousModeRecoveryParameterRow]] = {}
    for row in parameter_rows:
        grouped.setdefault((row.case_id, row.parameter), []).append(row)
    comparison_rows: list[ContinuousModeRecoveryParameterComparisonRow] = []
    for (_, parameter), rows in sorted(grouped.items()):
        by_engine = {row.recovery_engine: row for row in rows}
        if "bijux" not in by_engine or "geiger" not in by_engine:
            continue
        bijux_row = by_engine["bijux"]
        geiger_row = by_engine["geiger"]
        if abs(bijux_row.absolute_error - geiger_row.absolute_error) <= 1e-12:
            closer_engine = "tie"
        elif bijux_row.absolute_error < geiger_row.absolute_error:
            closer_engine = "bijux"
        else:
            closer_engine = "geiger"
        comparison_rows.append(
            ContinuousModeRecoveryParameterComparisonRow(
                case_id=bijux_row.case_id,
                generating_model=bijux_row.generating_model,
                parameter=parameter,
                true_value=bijux_row.true_value,
                bijux_estimated_value=bijux_row.estimated_value,
                geiger_estimated_value=geiger_row.estimated_value,
                bijux_absolute_error=bijux_row.absolute_error,
                geiger_absolute_error=geiger_row.absolute_error,
                bijux_within_tolerance=bijux_row.within_tolerance,
                geiger_within_tolerance=geiger_row.within_tolerance,
                closer_engine=closer_engine,
                tolerance=bijux_row.tolerance,
                interpretation=(
                    "Bijux and stored geiger recovery are compared directly against the same known simulation truth."
                ),
            )
        )
    return comparison_rows


def _parameter_row(
    *,
    case_id: str,
    generating_model: str,
    recovery_engine: str,
    fitted_model: str,
    fit_status: str,
    parameter: str,
    true_value: float,
    estimated_value: float,
    tolerance: float,
) -> ContinuousModeRecoveryParameterRow:
    absolute_error = abs(estimated_value - true_value)
    relative_error = 0.0 if true_value == 0.0 else absolute_error / abs(true_value)
    return ContinuousModeRecoveryParameterRow(
        case_id=case_id,
        generating_model=generating_model,
        recovery_engine=recovery_engine,
        fitted_model=fitted_model,
        fit_status=fit_status,
        parameter=parameter,
        true_value=true_value,
        estimated_value=estimated_value,
        absolute_error=absolute_error,
        relative_error=relative_error,
        tolerance=tolerance,
        within_tolerance=absolute_error <= tolerance,
        interpretation=(
            "The fitted parameter is compared directly against the known simulation truth for this case."
        ),
    )


def _truth_parameters(scenario: ContinuousModeRecoveryScenario) -> dict[str, float]:
    truths = {"sigma_squared": scenario.sigma**2}
    if scenario.generating_model == "ornstein-uhlenbeck":
        truths["alpha"] = scenario.alpha or 0.0
        if scenario.theta is not None:
            truths["theta"] = scenario.theta
    elif scenario.generating_model == "early-burst":
        truths["rate_change"] = scenario.rate_change or 0.0
    elif scenario.generating_model == "pagel-lambda":
        truths["lambda"] = scenario.lambda_value or 0.0
    elif scenario.generating_model == "pagel-kappa":
        truths["kappa"] = scenario.kappa or 0.0
    elif scenario.generating_model == "pagel-delta":
        truths["delta"] = scenario.delta or 0.0
    return truths


def _estimate_for_parameter(
    snapshot: _RecoveryFitSnapshot,
    parameter: str,
) -> float | None:
    if parameter == "sigma_squared":
        return snapshot.rate
    if parameter == "theta":
        return snapshot.root_state
    if snapshot.parameter_name == parameter:
        return snapshot.parameter_value
    return None


def _selection_match(
    expected_selected_model: str | None,
    selected_model: str | None,
) -> bool | None:
    if expected_selected_model is None:
        return None
    return expected_selected_model == selected_model
