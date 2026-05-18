from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from bijux_phylogenetics.comparative.evolutionary_modes import (
    ContinuousEvolutionaryModeFitReport,
    compare_fitcontinuous_model_ranking,
    fit_continuous_evolutionary_mode,
    transform_tree_for_evolutionary_mode,
)
from bijux_phylogenetics.comparative.geiger_fitcontinuous_recovery_reference import (
    GEIGER_FITCONTINUOUS_RECOVERY_REFERENCE_PAYLOADS,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.simulation import (
    ContinuousTraitSimulationReport,
    simulate_brownian_traits,
    simulate_early_burst_traits,
    simulate_ou_traits,
    write_continuous_trait_table,
)

_TRAIT_NAME = "value"
_LEGACY_CANDIDATE_MODES = ("brownian", "ornstein-uhlenbeck", "early-burst")
_FULL_CANDIDATE_MODES = (
    "brownian",
    "white-noise",
    "pagel-lambda",
    "pagel-kappa",
    "pagel-delta",
    "ornstein-uhlenbeck",
    "early-burst",
)
_GEIGER_MODEL_NAMES = {
    "brownian": "BM",
    "white-noise": "white",
    "pagel-lambda": "lambda",
    "pagel-kappa": "kappa",
    "pagel-delta": "delta",
    "ornstein-uhlenbeck": "OU",
    "early-burst": "EB",
}


@dataclass(slots=True)
class ContinuousModeRecoveryScenario:
    """One governed simulation-recovery case for a continuous `fitContinuous` lane."""

    case_id: str
    label: str
    generating_model: str
    expected_selected_model: str | None
    root_state: float
    sigma: float
    seed: int
    tree_path: Path | None = None
    alpha: float | None = None
    theta: float | None = None
    rate_change: float | None = None
    lambda_value: float | None = None
    kappa: float | None = None
    delta: float | None = None
    candidate_modes: tuple[str, ...] = _LEGACY_CANDIDATE_MODES
    lambda_bounds: tuple[float, float] = (0.0, 1.0)
    kappa_bounds: tuple[float, float] = (0.0, 3.0)
    delta_bounds: tuple[float, float] = (0.0, 3.0)
    ou_bounds: tuple[float, float] = (0.0, 10.0)
    early_burst_bounds: tuple[float, float] = (0.0, 50.0)
    parameter_tolerances: dict[str, float] = field(default_factory=dict)
    expected_warning_kinds: list[str] = field(default_factory=list)
    notes: str = ""


@dataclass(slots=True)
class ContinuousModeRecoveryParameterRow:
    """One truth-versus-fit parameter-recovery row for one recovery engine."""

    case_id: str
    generating_model: str
    recovery_engine: str
    fitted_model: str
    fit_status: str
    parameter: str
    true_value: float
    estimated_value: float
    absolute_error: float
    relative_error: float
    tolerance: float
    within_tolerance: bool
    interpretation: str


@dataclass(slots=True)
class ContinuousModeRecoveryParameterComparisonRow:
    """One paired Bijux-versus-geiger truth-recovery comparison for one parameter."""

    case_id: str
    generating_model: str
    parameter: str
    true_value: float
    bijux_estimated_value: float
    geiger_estimated_value: float
    bijux_absolute_error: float
    geiger_absolute_error: float
    bijux_within_tolerance: bool
    geiger_within_tolerance: bool
    closer_engine: str
    tolerance: float
    interpretation: str


@dataclass(slots=True)
class ContinuousModeRecoveryModelChoiceRow:
    """One candidate-model row from a governed model-comparison recovery review."""

    case_id: str
    generating_model: str
    recovery_engine: str
    expected_selected_model: str | None
    model: str
    parameter_count: int
    log_likelihood: float
    aic: float
    aicc: float
    selected: bool


@dataclass(slots=True)
class ContinuousModeRecoveryExecutionRow:
    """One execution-status row for either a fit or model-comparison review."""

    case_id: str
    recovery_engine: str
    operation: str
    fitted_model: str
    fit_status: str
    selected_model: str | None
    optimizer_name: str | None
    converged: bool | None
    hit_lower_parameter_boundary: bool | None
    hit_upper_parameter_boundary: bool | None
    warning_count: int
    failure_reason: str | None


@dataclass(slots=True)
class ContinuousModeRecoveryWarningRow:
    """One identifiability or weak-recovery warning observed during review."""

    case_id: str
    recovery_engine: str
    fitted_model: str
    kind: str
    message: str


@dataclass(slots=True)
class ContinuousModeRecoveryCaseReport:
    """Full recovery review for one deterministic simulation case."""

    scenario: ContinuousModeRecoveryScenario
    tree_path: Path
    traits_path: Path | None
    simulation: ContinuousTraitSimulationReport
    parameter_rows: list[ContinuousModeRecoveryParameterRow]
    parameter_comparison_rows: list[ContinuousModeRecoveryParameterComparisonRow]
    model_choice_rows: list[ContinuousModeRecoveryModelChoiceRow]
    execution_rows: list[ContinuousModeRecoveryExecutionRow]
    warning_rows: list[ContinuousModeRecoveryWarningRow]
    selected_model: str | None
    geiger_selected_model: str | None
    selection_matches_expectation: bool | None
    geiger_selection_matches_expectation: bool | None
    expected_warning_kinds_present: bool


@dataclass(slots=True)
class ContinuousModeRecoveryReport:
    """Integrated simulation-recovery benchmark over governed continuous cases."""

    default_tree_path: Path
    case_reports: list[ContinuousModeRecoveryCaseReport]


@dataclass(slots=True)
class _RecoveryFitSnapshot:
    engine: str
    fitted_model: str
    fit_status: str
    failure_reason: str | None
    parameter_name: str | None
    parameter_value: float | None
    rate: float | None
    optimizer_name: str | None
    converged: bool | None
    hit_lower_parameter_boundary: bool | None
    hit_upper_parameter_boundary: bool | None
    warning_rows: list[ContinuousModeRecoveryWarningRow]


@dataclass(slots=True)
class _RecoveryComparisonSnapshot:
    engine: str
    fit_status: str
    failure_reason: str | None
    selected_model: str | None
    rows: list[ContinuousModeRecoveryModelChoiceRow]
    warning_count: int


def run_continuous_mode_recovery(
    tree_path: Path,
    scenarios: list[ContinuousModeRecoveryScenario],
    *,
    artifacts_root: Path | None = None,
) -> ContinuousModeRecoveryReport:
    """Simulate, refit, and compare governed continuous-mode recovery cases."""
    if artifacts_root is not None:
        artifacts_root.mkdir(parents=True, exist_ok=True)
        case_reports = [
            _run_case(
                default_tree_path=tree_path,
                scenario=scenario,
                working_root=artifacts_root,
                persist_traits=True,
            )
            for scenario in scenarios
        ]
        return ContinuousModeRecoveryReport(
            default_tree_path=tree_path,
            case_reports=case_reports,
        )
    with TemporaryDirectory(prefix="continuous-mode-recovery-") as temporary_root:
        temporary_path = Path(temporary_root)
        case_reports = [
            _run_case(
                default_tree_path=tree_path,
                scenario=scenario,
                working_root=temporary_path,
                persist_traits=False,
            )
            for scenario in scenarios
        ]
    return ContinuousModeRecoveryReport(
        default_tree_path=tree_path,
        case_reports=case_reports,
    )


def write_continuous_mode_recovery_summary_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write one integrated recovery summary row per governed simulation case."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "label",
            "tree_path",
            "generating_model",
            "expected_selected_model",
            "bijux_selected_model",
            "geiger_selected_model",
            "bijux_selection_matches_expectation",
            "geiger_selection_matches_expectation",
            "parameter_row_count",
            "parameter_comparison_row_count",
            "bijux_parameter_pass_count",
            "geiger_parameter_pass_count",
            "parameter_closer_to_truth_count_bijux",
            "parameter_closer_to_truth_count_geiger",
            "expected_warning_count",
            "expected_warning_kinds_present",
            "warning_count",
            "notes",
        ],
        rows=[
            {
                "case_id": case.scenario.case_id,
                "label": case.scenario.label,
                "tree_path": str(case.tree_path),
                "generating_model": case.scenario.generating_model,
                "expected_selected_model": case.scenario.expected_selected_model or "",
                "bijux_selected_model": case.selected_model or "",
                "geiger_selected_model": case.geiger_selected_model or "",
                "bijux_selection_matches_expectation": _format_optional_bool(
                    case.selection_matches_expectation
                ),
                "geiger_selection_matches_expectation": _format_optional_bool(
                    case.geiger_selection_matches_expectation
                ),
                "parameter_row_count": str(len(case.parameter_rows)),
                "parameter_comparison_row_count": str(
                    len(case.parameter_comparison_rows)
                ),
                "bijux_parameter_pass_count": str(
                    sum(
                        1
                        for row in case.parameter_rows
                        if row.recovery_engine == "bijux" and row.within_tolerance
                    )
                ),
                "geiger_parameter_pass_count": str(
                    sum(
                        1
                        for row in case.parameter_rows
                        if row.recovery_engine == "geiger" and row.within_tolerance
                    )
                ),
                "parameter_closer_to_truth_count_bijux": str(
                    sum(
                        1
                        for row in case.parameter_comparison_rows
                        if row.closer_engine == "bijux"
                    )
                ),
                "parameter_closer_to_truth_count_geiger": str(
                    sum(
                        1
                        for row in case.parameter_comparison_rows
                        if row.closer_engine == "geiger"
                    )
                ),
                "expected_warning_count": str(
                    len(case.scenario.expected_warning_kinds)
                ),
                "expected_warning_kinds_present": str(
                    case.expected_warning_kinds_present
                ).lower(),
                "warning_count": str(len(case.warning_rows)),
                "notes": case.scenario.notes,
            }
            for case in report.case_reports
        ],
    )


def write_continuous_mode_recovery_parameter_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write all truth-versus-fit parameter-recovery rows across the cases."""
    rows = [
        {
            "case_id": row.case_id,
            "generating_model": row.generating_model,
            "recovery_engine": row.recovery_engine,
            "fitted_model": row.fitted_model,
            "fit_status": row.fit_status,
            "parameter": row.parameter,
            "true_value": _format_number(row.true_value),
            "estimated_value": _format_number(row.estimated_value),
            "absolute_error": _format_number(row.absolute_error),
            "relative_error": _format_number(row.relative_error),
            "tolerance": _format_number(row.tolerance),
            "within_tolerance": str(row.within_tolerance).lower(),
            "interpretation": row.interpretation,
        }
        for case in report.case_reports
        for row in case.parameter_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "generating_model",
            "recovery_engine",
            "fitted_model",
            "fit_status",
            "parameter",
            "true_value",
            "estimated_value",
            "absolute_error",
            "relative_error",
            "tolerance",
            "within_tolerance",
            "interpretation",
        ],
        rows=rows,
    )


def write_continuous_mode_recovery_parameter_comparison_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write paired Bijux-versus-geiger parameter-recovery comparisons."""
    rows = [
        {
            "case_id": row.case_id,
            "generating_model": row.generating_model,
            "parameter": row.parameter,
            "true_value": _format_number(row.true_value),
            "bijux_estimated_value": _format_number(row.bijux_estimated_value),
            "geiger_estimated_value": _format_number(row.geiger_estimated_value),
            "bijux_absolute_error": _format_number(row.bijux_absolute_error),
            "geiger_absolute_error": _format_number(row.geiger_absolute_error),
            "bijux_within_tolerance": str(row.bijux_within_tolerance).lower(),
            "geiger_within_tolerance": str(row.geiger_within_tolerance).lower(),
            "closer_engine": row.closer_engine,
            "tolerance": _format_number(row.tolerance),
            "interpretation": row.interpretation,
        }
        for case in report.case_reports
        for row in case.parameter_comparison_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "generating_model",
            "parameter",
            "true_value",
            "bijux_estimated_value",
            "geiger_estimated_value",
            "bijux_absolute_error",
            "geiger_absolute_error",
            "bijux_within_tolerance",
            "geiger_within_tolerance",
            "closer_engine",
            "tolerance",
            "interpretation",
        ],
        rows=rows,
    )


def write_continuous_mode_recovery_model_choice_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write candidate-model review rows for both recovery engines."""
    rows = [
        {
            "case_id": row.case_id,
            "generating_model": row.generating_model,
            "recovery_engine": row.recovery_engine,
            "expected_selected_model": row.expected_selected_model or "",
            "model": row.model,
            "parameter_count": str(row.parameter_count),
            "log_likelihood": _format_number(row.log_likelihood),
            "aic": _format_number(row.aic),
            "aicc": _format_number(row.aicc),
            "selected": str(row.selected).lower(),
        }
        for case in report.case_reports
        for row in case.model_choice_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "generating_model",
            "recovery_engine",
            "expected_selected_model",
            "model",
            "parameter_count",
            "log_likelihood",
            "aic",
            "aicc",
            "selected",
        ],
        rows=rows,
    )


def write_continuous_mode_recovery_execution_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write fit and model-comparison execution status rows for each engine."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "recovery_engine",
            "operation",
            "fitted_model",
            "fit_status",
            "selected_model",
            "optimizer_name",
            "converged",
            "hit_lower_parameter_boundary",
            "hit_upper_parameter_boundary",
            "warning_count",
            "failure_reason",
        ],
        rows=[
            {
                "case_id": row.case_id,
                "recovery_engine": row.recovery_engine,
                "operation": row.operation,
                "fitted_model": row.fitted_model,
                "fit_status": row.fit_status,
                "selected_model": row.selected_model or "",
                "optimizer_name": row.optimizer_name or "",
                "converged": _format_optional_bool(row.converged),
                "hit_lower_parameter_boundary": _format_optional_bool(
                    row.hit_lower_parameter_boundary
                ),
                "hit_upper_parameter_boundary": _format_optional_bool(
                    row.hit_upper_parameter_boundary
                ),
                "warning_count": str(row.warning_count),
                "failure_reason": row.failure_reason or "",
            }
            for case in report.case_reports
            for row in case.execution_rows
        ],
    )


def write_continuous_mode_recovery_warning_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write the governed warning rows observed across recovery cases."""
    return write_taxon_rows(
        path,
        columns=["case_id", "recovery_engine", "fitted_model", "kind", "message"],
        rows=[
            {
                "case_id": row.case_id,
                "recovery_engine": row.recovery_engine,
                "fitted_model": row.fitted_model,
                "kind": row.kind,
                "message": row.message,
            }
            for case in report.case_reports
            for row in case.warning_rows
        ],
    )


def _run_case(
    *,
    default_tree_path: Path,
    scenario: ContinuousModeRecoveryScenario,
    working_root: Path,
    persist_traits: bool,
) -> ContinuousModeRecoveryCaseReport:
    case_tree_path = default_tree_path if scenario.tree_path is None else scenario.tree_path
    simulation = _simulate_case(case_tree_path, scenario, working_root)
    traits_path = write_continuous_trait_table(
        working_root / f"{scenario.case_id}-traits.tsv",
        simulation,
    )
    bijux_fit = _build_bijux_fit_snapshot(case_tree_path, traits_path, scenario)
    bijux_comparison = _build_bijux_comparison_snapshot(
        case_tree_path,
        traits_path,
        scenario,
    )
    geiger_payload = GEIGER_FITCONTINUOUS_RECOVERY_REFERENCE_PAYLOADS[scenario.case_id]
    geiger_fit = _build_geiger_fit_snapshot(scenario, geiger_payload["fit_summary"])
    geiger_comparison = _build_geiger_comparison_snapshot(
        scenario,
        geiger_payload["comparison_summary"],
        geiger_payload["comparison_rows"],
    )
    parameter_rows = _build_parameter_rows(
        scenario=scenario,
        bijux_fit=bijux_fit,
        geiger_fit=geiger_fit,
    )
    parameter_comparison_rows = _build_parameter_comparison_rows(parameter_rows)
    warning_rows = [
        *bijux_fit.warning_rows,
        *geiger_fit.warning_rows,
    ]
    observed_warning_kinds = {
        row.kind
        for row in warning_rows
        if row.recovery_engine == "bijux"
    }
    execution_rows = [
        _build_execution_row(
            scenario.case_id,
            bijux_fit,
            operation="fit",
            selected_model=None,
        ),
        _build_execution_row(
            scenario.case_id,
            geiger_fit,
            operation="fit",
            selected_model=None,
        ),
        _build_comparison_execution_row(
            scenario.case_id,
            bijux_comparison,
        ),
        _build_comparison_execution_row(
            scenario.case_id,
            geiger_comparison,
        ),
    ]
    model_choice_rows = [
        *bijux_comparison.rows,
        *geiger_comparison.rows,
    ]
    return ContinuousModeRecoveryCaseReport(
        scenario=scenario,
        tree_path=case_tree_path,
        traits_path=traits_path if persist_traits else None,
        simulation=simulation,
        parameter_rows=parameter_rows,
        parameter_comparison_rows=parameter_comparison_rows,
        model_choice_rows=model_choice_rows,
        execution_rows=execution_rows,
        warning_rows=warning_rows,
        selected_model=bijux_comparison.selected_model,
        geiger_selected_model=geiger_comparison.selected_model,
        selection_matches_expectation=_selection_match(
            scenario.expected_selected_model,
            bijux_comparison.selected_model,
        ),
        geiger_selection_matches_expectation=_selection_match(
            scenario.expected_selected_model,
            geiger_comparison.selected_model,
        ),
        expected_warning_kinds_present=all(
            kind in observed_warning_kinds for kind in scenario.expected_warning_kinds
        ),
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
        trait=_TRAIT_NAME,
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


def _build_bijux_comparison_snapshot(
    tree_path: Path,
    traits_path: Path,
    scenario: ContinuousModeRecoveryScenario,
) -> _RecoveryComparisonSnapshot:
    report = compare_fitcontinuous_model_ranking(
        tree_path,
        traits_path,
        trait=_TRAIT_NAME,
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
    if snapshot.parameter_name == parameter:
        return snapshot.parameter_value
    return None


def _scenario_transform_parameter_value(scenario: ContinuousModeRecoveryScenario) -> float:
    if scenario.generating_model == "pagel-lambda":
        return scenario.lambda_value or 0.0
    if scenario.generating_model == "pagel-kappa":
        return scenario.kappa or 0.0
    if scenario.generating_model == "pagel-delta":
        return scenario.delta or 0.0
    raise ValueError(
        f"generating model does not use a branch-transform parameter: {scenario.generating_model}"
    )


def _selection_match(
    expected_selected_model: str | None,
    selected_model: str | None,
) -> bool | None:
    if expected_selected_model is None:
        return None
    return expected_selected_model == selected_model


def _format_number(value: float) -> str:
    return format(value, ".15g")


def _format_optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return str(value).lower()


def _optional_string(value: object) -> str | None:
    return value if isinstance(value, str) and value else None


def _optional_float(value: object) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _optional_bool(value: object) -> bool | None:
    if isinstance(value, bool):
        return value
    return None


def _geiger_candidate_model_names(
    candidate_modes: tuple[str, ...],
) -> list[str]:
    return [_GEIGER_MODEL_NAMES[mode] for mode in candidate_modes]


def geiger_fitcontinuous_recovery_reference_payload(case_id: str) -> dict[str, object]:
    """Expose the governed stored geiger recovery payload for one recovery case."""
    return GEIGER_FITCONTINUOUS_RECOVERY_REFERENCE_PAYLOADS[case_id]


def write_geiger_fitcontinuous_recovery_reference_payload_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write the stored governed geiger payload summaries used in the benchmark."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "fit_summary_json",
            "comparison_summary_json",
        ],
        rows=[
            {
                "case_id": case.scenario.case_id,
                "fit_summary_json": json.dumps(
                    GEIGER_FITCONTINUOUS_RECOVERY_REFERENCE_PAYLOADS[
                        case.scenario.case_id
                    ]["fit_summary"],
                    sort_keys=True,
                ),
                "comparison_summary_json": json.dumps(
                    GEIGER_FITCONTINUOUS_RECOVERY_REFERENCE_PAYLOADS[
                        case.scenario.case_id
                    ]["comparison_summary"],
                    sort_keys=True,
                ),
            }
            for case in report.case_reports
        ],
    )
