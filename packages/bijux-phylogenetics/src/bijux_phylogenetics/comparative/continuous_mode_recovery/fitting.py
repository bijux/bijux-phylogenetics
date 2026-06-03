from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.continuous_mode_recovery.reference_payloads import (
    GEIGER_FITCONTINUOUS_RECOVERY_REFERENCE_PAYLOADS,
)
from bijux_phylogenetics.comparative.evolutionary_modes import (
    ContinuousEvolutionaryModeFitReport,
    fit_continuous_evolutionary_mode,
    transform_tree_for_evolutionary_mode,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.simulation import (
    ContinuousTraitSimulationReport,
    simulate_brownian_traits,
    simulate_early_burst_traits,
    simulate_ou_traits,
)

from .formatting import _optional_bool, _optional_float, _optional_string
from .models import (
    RECOVERY_TRAIT_NAME,
    ContinuousModeRecoveryExecutionRow,
    ContinuousModeRecoveryScenario,
    ContinuousModeRecoveryWarningRow,
    _RecoveryComparisonSnapshot,
    _RecoveryFitSnapshot,
)


def _simulate_case(
    tree_path: Path,
    scenario: ContinuousModeRecoveryScenario,
    working_root: Path,
) -> ContinuousTraitSimulationReport:
    if scenario.generating_model == "brownian":
        return simulate_brownian_traits(
            tree_path,
            root_state=scenario.root_state,
            sigma=scenario.sigma,
            seed=scenario.seed,
        )
    if scenario.generating_model == "ornstein-uhlenbeck":
        return simulate_ou_traits(
            tree_path,
            root_state=scenario.root_state,
            sigma=scenario.sigma,
            alpha=scenario.alpha or 0.0,
            theta=scenario.theta or 0.0,
            seed=scenario.seed,
        )
    if scenario.generating_model == "early-burst":
        return simulate_early_burst_traits(
            tree_path,
            root_state=scenario.root_state,
            sigma=scenario.sigma,
            rate_change=scenario.rate_change or 0.0,
            seed=scenario.seed,
        )
    transformed_tree = transform_tree_for_evolutionary_mode(
        load_tree(tree_path),
        mode=scenario.generating_model,
        parameter_value=_scenario_transform_parameter_value(scenario),
    )
    transformed_tree_path = working_root / f"{scenario.case_id}-transformed-tree.nwk"
    transformed_tree_path.write_text(
        dumps_newick(transformed_tree) + "\n",
        encoding="utf-8",
    )
    return simulate_brownian_traits(
        transformed_tree_path,
        root_state=scenario.root_state,
        sigma=scenario.sigma,
        seed=scenario.seed,
    )


def _build_bijux_fit_snapshot(
    tree_path: Path,
    traits_path: Path,
    scenario: ContinuousModeRecoveryScenario,
) -> _RecoveryFitSnapshot:
    fit = fit_continuous_evolutionary_mode(
        tree_path,
        traits_path,
        trait=RECOVERY_TRAIT_NAME,
        mode=scenario.generating_model,
        lambda_bounds=scenario.lambda_bounds,
        kappa_bounds=scenario.kappa_bounds,
        delta_bounds=scenario.delta_bounds,
        ou_bounds=scenario.ou_bounds,
        early_burst_bounds=scenario.early_burst_bounds,
    )
    return _recovery_fit_snapshot_from_report(
        scenario.case_id,
        engine="bijux",
        fitted_model=scenario.generating_model,
        report=fit,
    )


def _recovery_fit_snapshot_from_report(
    case_id: str,
    *,
    engine: str,
    fitted_model: str,
    report: ContinuousEvolutionaryModeFitReport,
) -> _RecoveryFitSnapshot:
    warning_rows = [
        ContinuousModeRecoveryWarningRow(
            case_id=case_id,
            recovery_engine=engine,
            fitted_model=fitted_model,
            kind=warning.kind,
            message=warning.message,
        )
        for warning in report.identifiability_warnings
    ]
    return _RecoveryFitSnapshot(
        engine=engine,
        fitted_model=fitted_model,
        fit_status="ok",
        failure_reason=None,
        parameter_name=report.parameter_name,
        parameter_value=report.parameter_value,
        root_state=report.root_state,
        rate=report.rate,
        optimizer_name=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.optimizer_name
        ),
        converged=(
            True
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.converged
        ),
        hit_lower_parameter_boundary=(
            False
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.hit_lower_boundary
        ),
        hit_upper_parameter_boundary=(
            False
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.hit_upper_boundary
        ),
        warning_rows=warning_rows,
    )


def _build_geiger_fit_snapshot(
    scenario: ContinuousModeRecoveryScenario,
    summary: dict[str, object],
) -> _RecoveryFitSnapshot:
    warning_rows: list[ContinuousModeRecoveryWarningRow] = []
    if bool(summary.get("hit_lower_parameter_boundary")):
        warning_rows.append(
            ContinuousModeRecoveryWarningRow(
                case_id=scenario.case_id,
                recovery_engine="geiger",
                fitted_model=scenario.generating_model,
                kind="boundary_lower",
                message="The stored geiger fit hit the lower governed parameter boundary.",
            )
        )
    if bool(summary.get("hit_upper_parameter_boundary")):
        warning_rows.append(
            ContinuousModeRecoveryWarningRow(
                case_id=scenario.case_id,
                recovery_engine="geiger",
                fitted_model=scenario.generating_model,
                kind="boundary_upper",
                message="The stored geiger fit hit the upper governed parameter boundary.",
            )
        )
    optimizer_result = summary.get("optimizer_result")
    converged = None
    optimizer_name = None
    if isinstance(optimizer_result, dict):
        optimizer_name = _optional_string(optimizer_result.get("best_method"))
        convergence_code = optimizer_result.get("convergence_code")
        if isinstance(convergence_code, int):
            converged = convergence_code == 0
    return _RecoveryFitSnapshot(
        engine="geiger",
        fitted_model=scenario.generating_model,
        fit_status="ok",
        failure_reason=None,
        parameter_name=_optional_string(summary.get("parameter_name")),
        parameter_value=_optional_float(summary.get("parameter_value")),
        root_state=_optional_float(summary.get("root_state")),
        rate=_optional_float(summary.get("rate")),
        optimizer_name=optimizer_name,
        converged=converged,
        hit_lower_parameter_boundary=_optional_bool(
            summary.get("hit_lower_parameter_boundary")
        ),
        hit_upper_parameter_boundary=_optional_bool(
            summary.get("hit_upper_parameter_boundary")
        ),
        warning_rows=warning_rows,
    )


def _build_execution_row(
    case_id: str,
    snapshot: _RecoveryFitSnapshot,
    *,
    operation: str,
    selected_model: str | None,
) -> ContinuousModeRecoveryExecutionRow:
    return ContinuousModeRecoveryExecutionRow(
        case_id=case_id,
        recovery_engine=snapshot.engine,
        operation=operation,
        fitted_model=snapshot.fitted_model,
        fit_status=snapshot.fit_status,
        selected_model=selected_model,
        optimizer_name=snapshot.optimizer_name,
        converged=snapshot.converged,
        hit_lower_parameter_boundary=snapshot.hit_lower_parameter_boundary,
        hit_upper_parameter_boundary=snapshot.hit_upper_parameter_boundary,
        warning_count=len(snapshot.warning_rows),
        failure_reason=snapshot.failure_reason,
    )


def _build_comparison_execution_row(
    case_id: str,
    snapshot: _RecoveryComparisonSnapshot,
) -> ContinuousModeRecoveryExecutionRow:
    return ContinuousModeRecoveryExecutionRow(
        case_id=case_id,
        recovery_engine=snapshot.engine,
        operation="model-comparison",
        fitted_model="candidate-set",
        fit_status=snapshot.fit_status,
        selected_model=snapshot.selected_model,
        optimizer_name=None,
        converged=None,
        hit_lower_parameter_boundary=None,
        hit_upper_parameter_boundary=None,
        warning_count=snapshot.warning_count,
        failure_reason=snapshot.failure_reason,
    )


def _geiger_payload_for_case(case_id: str) -> dict[str, object]:
    return GEIGER_FITCONTINUOUS_RECOVERY_REFERENCE_PAYLOADS[case_id]


def _scenario_transform_parameter_value(
    scenario: ContinuousModeRecoveryScenario,
) -> float:
    if scenario.generating_model == "pagel-lambda":
        return scenario.lambda_value or 0.0
    if scenario.generating_model == "pagel-kappa":
        return scenario.kappa or 0.0
    if scenario.generating_model == "pagel-delta":
        return scenario.delta or 0.0
    raise ValueError(
        f"generating model does not use a branch-transform parameter: {scenario.generating_model}"
    )
