from __future__ import annotations

from .models import (
    DiscreteModeRecoveryModelChoiceRow,
    DiscreteModeRecoveryParameterComparisonRow,
    DiscreteModeRecoveryParameterRow,
    DiscreteModeRecoveryRateComparisonRow,
    DiscreteModeRecoveryRateRow,
    DiscreteModeRecoveryScenario,
    _DiscreteRecoveryFitSnapshot,
)


def _build_parameter_rows(
    *,
    scenario: DiscreteModeRecoveryScenario,
    bijux_fit: _DiscreteRecoveryFitSnapshot,
    geiger_fit: _DiscreteRecoveryFitSnapshot,
) -> list[DiscreteModeRecoveryParameterRow]:
    transform = scenario.transform or ""
    true_value = scenario.transform_parameter_value
    if not transform or true_value is None:
        return []
    rows: list[DiscreteModeRecoveryParameterRow] = []
    for engine, snapshot in (("bijux", bijux_fit), ("geiger", geiger_fit)):
        parameter = snapshot.parameter_name
        estimate = snapshot.parameter_value
        if parameter is None or estimate is None:
            continue
        tolerance = scenario.parameter_tolerances.get(parameter)
        if tolerance is None:
            continue
        rows.append(
            _parameter_row(
                case_id=scenario.case_id,
                generating_model=scenario.generating_model,
                transform=transform,
                recovery_engine=engine,
                fitted_model=snapshot.fitted_model,
                fit_status=snapshot.fit_status,
                parameter=parameter,
                true_value=true_value,
                estimated_value=estimate,
                tolerance=tolerance,
            )
        )
    return rows


def _build_parameter_comparison_rows(
    parameter_rows: list[DiscreteModeRecoveryParameterRow],
) -> list[DiscreteModeRecoveryParameterComparisonRow]:
    grouped: dict[tuple[str, str], list[DiscreteModeRecoveryParameterRow]] = {}
    for row in parameter_rows:
        grouped.setdefault((row.case_id, row.parameter), []).append(row)
    comparison_rows: list[DiscreteModeRecoveryParameterComparisonRow] = []
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
            DiscreteModeRecoveryParameterComparisonRow(
                case_id=bijux_row.case_id,
                generating_model=bijux_row.generating_model,
                transform=bijux_row.transform,
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
                    "Bijux and stored geiger transform recovery are compared directly against the same known simulation truth."
                ),
            )
        )
    return comparison_rows


def _parameter_row(
    *,
    case_id: str,
    generating_model: str,
    transform: str,
    recovery_engine: str,
    fitted_model: str,
    fit_status: str,
    parameter: str,
    true_value: float,
    estimated_value: float,
    tolerance: float,
) -> DiscreteModeRecoveryParameterRow:
    absolute_error = abs(estimated_value - true_value)
    relative_error = 0.0 if true_value == 0.0 else absolute_error / abs(true_value)
    return DiscreteModeRecoveryParameterRow(
        case_id=case_id,
        generating_model=generating_model,
        transform=transform,
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
            "The fitted transform parameter is compared directly against the known simulation truth for this discrete Mk case."
        ),
    )


def _build_rate_rows(
    *,
    scenario: DiscreteModeRecoveryScenario,
    bijux_fit: _DiscreteRecoveryFitSnapshot,
    geiger_fit: _DiscreteRecoveryFitSnapshot,
) -> list[DiscreteModeRecoveryRateRow]:
    tolerance = scenario.rate_tolerance
    truth_by_pair = {
        (row.source_state, row.target_state): row.rate for row in scenario.rate_rows
    }
    rows: list[DiscreteModeRecoveryRateRow] = []
    for engine, snapshot in (("bijux", bijux_fit), ("geiger", geiger_fit)):
        estimate_by_pair = {
            (row.source_state, row.target_state): row.rate for row in snapshot.rate_rows
        }
        for pair, true_rate in sorted(truth_by_pair.items()):
            estimated_rate = estimate_by_pair.get(pair)
            if estimated_rate is None:
                absolute_error = None
                relative_error = None
                within_tolerance = None
                interpretation = "The fitted transition rate is missing because one or more states were not observed after simulation pruning, so this review row records the omitted estimate explicitly."
            else:
                absolute_error = abs(estimated_rate - true_rate)
                relative_error = (
                    0.0 if true_rate == 0.0 else absolute_error / abs(true_rate)
                )
                within_tolerance = (
                    None if tolerance is None else absolute_error <= tolerance
                )
                interpretation = (
                    "The fitted transition rate is compared directly against the known simulation truth for this discrete Mk case."
                    if tolerance is not None
                    else "The fitted transition rate is compared directly against the known simulation truth for this discrete Mk review case without a governed pass-fail tolerance."
                )
            rows.append(
                DiscreteModeRecoveryRateRow(
                    case_id=scenario.case_id,
                    generating_model=scenario.generating_model,
                    transform=scenario.transform or "",
                    recovery_engine=engine,
                    fitted_model=snapshot.fitted_model,
                    fit_status=snapshot.fit_status,
                    source_state=pair[0],
                    target_state=pair[1],
                    true_rate=true_rate,
                    estimated_rate=estimated_rate,
                    absolute_error=absolute_error,
                    relative_error=relative_error,
                    tolerance=tolerance,
                    within_tolerance=within_tolerance,
                    interpretation=interpretation,
                )
            )
    return rows


def _build_rate_comparison_rows(
    rate_rows: list[DiscreteModeRecoveryRateRow],
) -> list[DiscreteModeRecoveryRateComparisonRow]:
    grouped: dict[tuple[str, str, str], list[DiscreteModeRecoveryRateRow]] = {}
    for row in rate_rows:
        grouped.setdefault(
            (row.case_id, row.source_state, row.target_state),
            [],
        ).append(row)
    comparison_rows: list[DiscreteModeRecoveryRateComparisonRow] = []
    for (_, source_state, target_state), rows in sorted(grouped.items()):
        by_engine = {row.recovery_engine: row for row in rows}
        if "bijux" not in by_engine or "geiger" not in by_engine:
            continue
        bijux_row = by_engine["bijux"]
        geiger_row = by_engine["geiger"]
        if (
            bijux_row.absolute_error is None
            or geiger_row.absolute_error is None
            or bijux_row.estimated_rate is None
            or geiger_row.estimated_rate is None
        ):
            closer_engine = "missing-estimate"
        elif abs(bijux_row.absolute_error - geiger_row.absolute_error) <= 1e-12:
            closer_engine = "tie"
        elif bijux_row.absolute_error < geiger_row.absolute_error:
            closer_engine = "bijux"
        else:
            closer_engine = "geiger"
        comparison_rows.append(
            DiscreteModeRecoveryRateComparisonRow(
                case_id=bijux_row.case_id,
                generating_model=bijux_row.generating_model,
                transform=bijux_row.transform,
                source_state=source_state,
                target_state=target_state,
                true_rate=bijux_row.true_rate,
                bijux_estimated_rate=bijux_row.estimated_rate,
                geiger_estimated_rate=geiger_row.estimated_rate,
                bijux_absolute_error=bijux_row.absolute_error,
                geiger_absolute_error=geiger_row.absolute_error,
                bijux_within_tolerance=bijux_row.within_tolerance,
                geiger_within_tolerance=geiger_row.within_tolerance,
                closer_engine=closer_engine,
                tolerance=bijux_row.tolerance,
                interpretation=(
                    "Bijux and stored geiger discrete recovery are compared directly against the same known transition-rate truth."
                    if closer_engine != "missing-estimate"
                    else "Bijux and stored geiger discrete recovery both preserve the omitted transition explicitly because at least one fitted surface did not estimate that state pair."
                ),
            )
        )
    return comparison_rows


def _selection_match(
    expected_selected_model: str | None,
    selected_model: str | None,
) -> bool | None:
    if expected_selected_model is None:
        return None
    return expected_selected_model == selected_model


def _owner_model_name(model_name: str) -> str:
    mapping = {
        "ER": "equal-rates",
        "SYM": "symmetric",
        "ARD": "all-rates-different",
    }
    if model_name in mapping:
        return mapping[model_name]
    return model_name


def _select_best_model(rows: list[DiscreteModeRecoveryModelChoiceRow]) -> str:
    ranked_rows = sorted(rows, key=lambda row: (row.aicc, row.aic, row.model))
    return ranked_rows[0].model
