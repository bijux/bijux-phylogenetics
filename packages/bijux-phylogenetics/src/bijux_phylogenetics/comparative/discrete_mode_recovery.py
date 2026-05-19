from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from bijux_phylogenetics.comparative.discrete_mk import (
    DiscreteMkFitReport,
    fit_discrete_mk_model,
)
from bijux_phylogenetics.comparative.geiger_fitdiscrete_recovery_reference import (
    GEIGER_FITDISCRETE_RECOVERY_REFERENCE_PAYLOADS,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.simulation import (
    DiscreteHistoryRateRow,
    DiscreteTraitSimulationReport,
    simulate_discrete_histories,
    write_discrete_trait_table,
)

_TRAIT_NAME = "state"
_CANDIDATE_MODELS = ("equal-rates", "symmetric", "all-rates-different")
_GEIGER_MODEL_NAMES = {
    "equal-rates": "ER",
    "symmetric": "SYM",
    "all-rates-different": "ARD",
}


@dataclass(slots=True)
class DiscreteModeRecoveryScenario:
    """One governed simulation-recovery case for discrete Mk fitting."""

    case_id: str
    label: str
    generating_model: str
    expected_selected_model: str | None
    states: list[str]
    rate_rows: list[DiscreteHistoryRateRow]
    root_state: str
    seed: int
    tree_path: Path | None = None
    candidate_models: tuple[str, ...] = _CANDIDATE_MODELS
    rate_tolerance: float | None = None
    expected_overparameterized: bool = False
    expected_warning_kinds: list[str] = field(default_factory=list)
    notes: str = ""


@dataclass(slots=True)
class DiscreteModeRecoveryRateRow:
    """One truth-versus-fit transition-rate recovery row for one engine."""

    case_id: str
    generating_model: str
    recovery_engine: str
    fitted_model: str
    fit_status: str
    source_state: str
    target_state: str
    true_rate: float
    estimated_rate: float | None
    absolute_error: float | None
    relative_error: float | None
    tolerance: float | None
    within_tolerance: bool | None
    interpretation: str


@dataclass(slots=True)
class DiscreteModeRecoveryRateComparisonRow:
    """One paired Bijux-versus-geiger transition-rate recovery comparison."""

    case_id: str
    generating_model: str
    source_state: str
    target_state: str
    true_rate: float
    bijux_estimated_rate: float | None
    geiger_estimated_rate: float | None
    bijux_absolute_error: float | None
    geiger_absolute_error: float | None
    bijux_within_tolerance: bool | None
    geiger_within_tolerance: bool | None
    closer_engine: str
    tolerance: float | None
    interpretation: str


@dataclass(slots=True)
class DiscreteModeRecoveryModelChoiceRow:
    """One candidate-model row from a discrete recovery model-comparison review."""

    case_id: str
    generating_model: str
    recovery_engine: str
    expected_selected_model: str | None
    model: str
    parameter_count: int
    log_likelihood: float
    aic: float
    aicc: float
    overparameterized: bool
    selected: bool


@dataclass(slots=True)
class DiscreteModeRecoveryExecutionRow:
    """One execution-status row for either a fit or model-comparison review."""

    case_id: str
    recovery_engine: str
    operation: str
    fitted_model: str
    fit_status: str
    selected_model: str | None
    optimizer_name: str | None
    converged: bool | None
    hit_lower_parameter_bound: bool | None
    hit_upper_parameter_bound: bool | None
    overparameterized: bool
    warning_count: int
    failure_reason: str | None


@dataclass(slots=True)
class DiscreteModeRecoveryWarningRow:
    """One weak-identifiability or fit-risk warning observed during review."""

    case_id: str
    recovery_engine: str
    fitted_model: str
    kind: str
    message: str


@dataclass(slots=True)
class DiscreteModeRecoveryCaseReport:
    """Full discrete simulation-recovery review for one governed case."""

    scenario: DiscreteModeRecoveryScenario
    tree_path: Path
    traits_path: Path | None
    simulation: DiscreteTraitSimulationReport
    rate_rows: list[DiscreteModeRecoveryRateRow]
    rate_comparison_rows: list[DiscreteModeRecoveryRateComparisonRow]
    model_choice_rows: list[DiscreteModeRecoveryModelChoiceRow]
    execution_rows: list[DiscreteModeRecoveryExecutionRow]
    warning_rows: list[DiscreteModeRecoveryWarningRow]
    selected_model: str | None
    geiger_selected_model: str | None
    selection_matches_expectation: bool | None
    geiger_selection_matches_expectation: bool | None
    overparameterized_review_matches_expectation: bool
    expected_warning_kinds_present: bool


@dataclass(slots=True)
class DiscreteModeRecoveryReport:
    """Integrated discrete Mk simulation-recovery benchmark report."""

    default_tree_path: Path
    case_reports: list[DiscreteModeRecoveryCaseReport]


@dataclass(slots=True)
class _DiscreteRecoveryFitSnapshot:
    engine: str
    fitted_model: str
    fit_status: str
    failure_reason: str | None
    selected_model: str | None
    parameter_count: int
    log_likelihood: float
    aic: float
    aicc: float
    overparameterized: bool
    optimizer_name: str | None
    converged: bool | None
    hit_lower_parameter_bound: bool | None
    hit_upper_parameter_bound: bool | None
    rate_rows: list[DiscreteHistoryRateRow]
    warning_rows: list[DiscreteModeRecoveryWarningRow]


def run_discrete_mode_recovery(
    tree_path: Path,
    scenarios: list[DiscreteModeRecoveryScenario],
    *,
    artifacts_root: Path | None = None,
) -> DiscreteModeRecoveryReport:
    """Simulate, refit, and compare governed discrete Mk recovery cases."""
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
        return DiscreteModeRecoveryReport(
            default_tree_path=tree_path,
            case_reports=case_reports,
        )
    with TemporaryDirectory(prefix="discrete-mode-recovery-") as temporary_root:
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
    return DiscreteModeRecoveryReport(
        default_tree_path=tree_path,
        case_reports=case_reports,
    )


def write_discrete_mode_recovery_summary_table(
    path: Path,
    report: DiscreteModeRecoveryReport,
) -> Path:
    """Write one integrated summary row per discrete recovery case."""
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
            "overparameterized_review_matches_expectation",
            "rate_row_count",
            "rate_comparison_row_count",
            "bijux_rate_pass_count",
            "geiger_rate_pass_count",
            "rate_closer_to_truth_count_bijux",
            "rate_closer_to_truth_count_geiger",
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
                "overparameterized_review_matches_expectation": str(
                    case.overparameterized_review_matches_expectation
                ).lower(),
                "rate_row_count": str(len(case.rate_rows)),
                "rate_comparison_row_count": str(len(case.rate_comparison_rows)),
                "bijux_rate_pass_count": str(
                    sum(
                        1
                        for row in case.rate_rows
                        if row.recovery_engine == "bijux"
                        and row.within_tolerance is True
                    )
                ),
                "geiger_rate_pass_count": str(
                    sum(
                        1
                        for row in case.rate_rows
                        if row.recovery_engine == "geiger"
                        and row.within_tolerance is True
                    )
                ),
                "rate_closer_to_truth_count_bijux": str(
                    sum(
                        1
                        for row in case.rate_comparison_rows
                        if row.closer_engine == "bijux"
                    )
                ),
                "rate_closer_to_truth_count_geiger": str(
                    sum(
                        1
                        for row in case.rate_comparison_rows
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


def write_discrete_mode_recovery_rate_table(
    path: Path,
    report: DiscreteModeRecoveryReport,
) -> Path:
    """Write all transition-rate truth-versus-fit rows across the cases."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "generating_model",
            "recovery_engine",
            "fitted_model",
            "fit_status",
            "source_state",
            "target_state",
            "true_rate",
            "estimated_rate",
            "absolute_error",
            "relative_error",
            "tolerance",
            "within_tolerance",
            "interpretation",
        ],
        rows=[
            {
                "case_id": row.case_id,
                "generating_model": row.generating_model,
                "recovery_engine": row.recovery_engine,
                "fitted_model": row.fitted_model,
                "fit_status": row.fit_status,
                "source_state": row.source_state,
                "target_state": row.target_state,
                "true_rate": _format_number(row.true_rate),
                "estimated_rate": _format_optional_number(row.estimated_rate),
                "absolute_error": _format_optional_number(row.absolute_error),
                "relative_error": _format_optional_number(row.relative_error),
                "tolerance": _format_optional_number(row.tolerance),
                "within_tolerance": _format_optional_bool(row.within_tolerance),
                "interpretation": row.interpretation,
            }
            for case in report.case_reports
            for row in case.rate_rows
        ],
    )


def write_discrete_mode_recovery_rate_comparison_table(
    path: Path,
    report: DiscreteModeRecoveryReport,
) -> Path:
    """Write paired Bijux-versus-geiger transition-rate recovery comparisons."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "generating_model",
            "source_state",
            "target_state",
            "true_rate",
            "bijux_estimated_rate",
            "geiger_estimated_rate",
            "bijux_absolute_error",
            "geiger_absolute_error",
            "bijux_within_tolerance",
            "geiger_within_tolerance",
            "closer_engine",
            "tolerance",
            "interpretation",
        ],
        rows=[
            {
                "case_id": row.case_id,
                "generating_model": row.generating_model,
                "source_state": row.source_state,
                "target_state": row.target_state,
                "true_rate": _format_number(row.true_rate),
                "bijux_estimated_rate": _format_optional_number(
                    row.bijux_estimated_rate
                ),
                "geiger_estimated_rate": _format_optional_number(
                    row.geiger_estimated_rate
                ),
                "bijux_absolute_error": _format_optional_number(
                    row.bijux_absolute_error
                ),
                "geiger_absolute_error": _format_optional_number(
                    row.geiger_absolute_error
                ),
                "bijux_within_tolerance": _format_optional_bool(
                    row.bijux_within_tolerance
                ),
                "geiger_within_tolerance": _format_optional_bool(
                    row.geiger_within_tolerance
                ),
                "closer_engine": row.closer_engine,
                "tolerance": _format_optional_number(row.tolerance),
                "interpretation": row.interpretation,
            }
            for case in report.case_reports
            for row in case.rate_comparison_rows
        ],
    )


def write_discrete_mode_recovery_model_choice_table(
    path: Path,
    report: DiscreteModeRecoveryReport,
) -> Path:
    """Write candidate-model selection review rows for both engines."""
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
            "overparameterized",
            "selected",
        ],
        rows=[
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
                "overparameterized": str(row.overparameterized).lower(),
                "selected": str(row.selected).lower(),
            }
            for case in report.case_reports
            for row in case.model_choice_rows
        ],
    )


def write_discrete_mode_recovery_execution_table(
    path: Path,
    report: DiscreteModeRecoveryReport,
) -> Path:
    """Write fit and model-comparison execution rows for each engine."""
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
            "hit_lower_parameter_bound",
            "hit_upper_parameter_bound",
            "overparameterized",
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
                "hit_lower_parameter_bound": _format_optional_bool(
                    row.hit_lower_parameter_bound
                ),
                "hit_upper_parameter_bound": _format_optional_bool(
                    row.hit_upper_parameter_bound
                ),
                "overparameterized": str(row.overparameterized).lower(),
                "warning_count": str(row.warning_count),
                "failure_reason": row.failure_reason or "",
            }
            for case in report.case_reports
            for row in case.execution_rows
        ],
    )


def write_discrete_mode_recovery_warning_table(
    path: Path,
    report: DiscreteModeRecoveryReport,
) -> Path:
    """Write all benchmark warning rows across the cases."""
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


def write_geiger_fitdiscrete_recovery_reference_payload_table(
    path: Path,
    report: DiscreteModeRecoveryReport,
) -> Path:
    """Write the stored governed geiger payload summaries used in the benchmark."""
    return write_taxon_rows(
        path,
        columns=["case_id", "fit_summary_json", "comparison_summary_json"],
        rows=[
            {
                "case_id": case.scenario.case_id,
                "fit_summary_json": json.dumps(
                    GEIGER_FITDISCRETE_RECOVERY_REFERENCE_PAYLOADS[
                        case.scenario.case_id
                    ]["fit_summary"],
                    sort_keys=True,
                ),
                "comparison_summary_json": json.dumps(
                    GEIGER_FITDISCRETE_RECOVERY_REFERENCE_PAYLOADS[
                        case.scenario.case_id
                    ]["comparison_summary"],
                    sort_keys=True,
                ),
            }
            for case in report.case_reports
        ],
    )


def geiger_fitdiscrete_recovery_reference_payload(case_id: str) -> dict[str, object]:
    """Expose the governed stored geiger recovery payload for one discrete case."""
    return GEIGER_FITDISCRETE_RECOVERY_REFERENCE_PAYLOADS[case_id]


def _run_case(
    *,
    default_tree_path: Path,
    scenario: DiscreteModeRecoveryScenario,
    working_root: Path,
    persist_traits: bool,
) -> DiscreteModeRecoveryCaseReport:
    case_tree_path = default_tree_path if scenario.tree_path is None else scenario.tree_path
    simulation = _simulate_case(case_tree_path, scenario)
    traits_path = write_discrete_trait_table(
        working_root / f"{scenario.case_id}-traits.tsv",
        simulation,
    )
    bijux_fit = _build_bijux_fit_snapshot(case_tree_path, traits_path, scenario)
    bijux_model_rows, bijux_selected_model = _build_bijux_model_choice_rows(
        case_tree_path,
        traits_path,
        scenario,
    )
    geiger_payload = GEIGER_FITDISCRETE_RECOVERY_REFERENCE_PAYLOADS[scenario.case_id]
    geiger_fit = _build_geiger_fit_snapshot(scenario, geiger_payload["fit_summary"], geiger_payload["fit_rows"])
    geiger_model_rows, geiger_selected_model = _build_geiger_model_choice_rows(
        scenario,
        geiger_payload["comparison_summary"],
        geiger_payload["comparison_rows"],
    )
    rate_rows = _build_rate_rows(scenario=scenario, bijux_fit=bijux_fit, geiger_fit=geiger_fit)
    rate_comparison_rows = _build_rate_comparison_rows(rate_rows)
    warning_rows = [*bijux_fit.warning_rows, *geiger_fit.warning_rows]
    observed_warning_kinds = {
        row.kind for row in warning_rows if row.recovery_engine == "bijux"
    }
    overparameterized_review_matches_expectation = (
        bijux_fit.overparameterized == scenario.expected_overparameterized
        and geiger_fit.overparameterized == scenario.expected_overparameterized
    )
    execution_rows = [
        _build_execution_row(scenario.case_id, bijux_fit, operation="fit"),
        _build_execution_row(scenario.case_id, geiger_fit, operation="fit"),
        _build_model_comparison_execution_row(
            scenario.case_id,
            "bijux",
            bijux_selected_model,
            bijux_model_rows,
        ),
        _build_model_comparison_execution_row(
            scenario.case_id,
            "geiger",
            geiger_selected_model,
            geiger_model_rows,
        ),
    ]
    return DiscreteModeRecoveryCaseReport(
        scenario=scenario,
        tree_path=case_tree_path,
        traits_path=traits_path if persist_traits else None,
        simulation=simulation,
        rate_rows=rate_rows,
        rate_comparison_rows=rate_comparison_rows,
        model_choice_rows=[*bijux_model_rows, *geiger_model_rows],
        execution_rows=execution_rows,
        warning_rows=warning_rows,
        selected_model=bijux_selected_model,
        geiger_selected_model=geiger_selected_model,
        selection_matches_expectation=_selection_match(
            scenario.expected_selected_model,
            bijux_selected_model,
        ),
        geiger_selection_matches_expectation=_selection_match(
            scenario.expected_selected_model,
            geiger_selected_model,
        ),
        overparameterized_review_matches_expectation=overparameterized_review_matches_expectation,
        expected_warning_kinds_present=all(
            kind in observed_warning_kinds for kind in scenario.expected_warning_kinds
        ),
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
        trait=_TRAIT_NAME,
        model=scenario.generating_model,
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
                message=(
                    "One or more fitted discrete Mk rates hit an optimizer bound."
                ),
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
        parameter_count=report.parameter_count,
        log_likelihood=report.log_likelihood,
        aic=report.aic,
        aicc=report.aicc,
        overparameterized=report.overparameterized,
        optimizer_name=report.optimizer_diagnostics.optimizer_name,
        converged=report.optimizer_diagnostics.converged,
        hit_lower_parameter_bound=report.optimizer_diagnostics.hit_lower_parameter_bound,
        hit_upper_parameter_bound=report.optimizer_diagnostics.hit_upper_parameter_bound,
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
    rows: list[DiscreteModeRecoveryModelChoiceRow] = []
    for model in scenario.candidate_models:
        report = fit_discrete_mk_model(
            tree_path,
            traits_path,
            trait=_TRAIT_NAME,
            model=model,
        )
        rows.append(
            DiscreteModeRecoveryModelChoiceRow(
                case_id=scenario.case_id,
                generating_model=scenario.generating_model,
                recovery_engine="bijux",
                expected_selected_model=scenario.expected_selected_model,
                model=model,
                parameter_count=report.parameter_count,
                log_likelihood=report.log_likelihood,
                aic=report.aic,
                aicc=report.aicc,
                overparameterized=report.overparameterized,
                selected=False,
            )
        )
    selected_model = _select_best_model(rows)
    for row in rows:
        row.selected = row.model == selected_model
    return rows, selected_model


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
        parameter_count=parameter_count,
        log_likelihood=float(summary["log_likelihood"]),
        aic=float(summary["aic"]),
        aicc=float(summary["aicc"]),
        overparameterized=parameter_count >= taxon_count,
        optimizer_name=optimizer_name,
        converged=converged,
        hit_lower_parameter_bound=None,
        hit_upper_parameter_bound=None,
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
            recovery_engine="geiger",
            expected_selected_model=scenario.expected_selected_model,
            model=_owner_model_name(str(row["model"])),
            parameter_count=int(row["parameter_count"]),
            log_likelihood=float(row["log_likelihood"]),
            aic=float(row["aic"]),
            aicc=float(row["aicc"]),
            overparameterized=int(row["parameter_count"]) >= int(summary["taxon_count"]),
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


def _build_rate_rows(
    *,
    scenario: DiscreteModeRecoveryScenario,
    bijux_fit: _DiscreteRecoveryFitSnapshot,
    geiger_fit: _DiscreteRecoveryFitSnapshot,
) -> list[DiscreteModeRecoveryRateRow]:
    tolerance = scenario.rate_tolerance
    truth_by_pair = {(row.source_state, row.target_state): row.rate for row in scenario.rate_rows}
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
                interpretation = (
                    "The fitted transition rate is missing because one or more states were not observed after simulation pruning, so this review row records the omitted estimate explicitly."
                )
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
                    recovery_engine=engine,
                    fitted_model=scenario.generating_model,
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
        grouped.setdefault((row.case_id, row.source_state, row.target_state), []).append(row)
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


def _select_best_model(rows: list[DiscreteModeRecoveryModelChoiceRow]) -> str:
    ranked_rows = sorted(rows, key=lambda row: (row.aicc, row.aic, row.model))
    return ranked_rows[0].model


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


def _format_number(value: float) -> str:
    return format(value, ".15g")


def _format_optional_number(value: float | None) -> str:
    if value is None:
        return ""
    return _format_number(value)


def _format_optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return str(value).lower()


def _optional_string(value: object) -> str | None:
    return value if isinstance(value, str) and value else None
