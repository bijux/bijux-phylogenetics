from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from tempfile import TemporaryDirectory

from bijux_phylogenetics.comparative.brownian_trait_evolution import (
    BrownianTraitEvolutionSummaryReport,
    summarize_brownian_trait_evolution,
)
from bijux_phylogenetics.comparative.early_burst_trait_evolution import (
    EarlyBurstTraitEvolutionSummaryReport,
    summarize_early_burst_trait_evolution,
)
from bijux_phylogenetics.comparative.evolutionary_modes import (
    ContinuousEvolutionaryModeComparisonReport,
    compare_continuous_evolutionary_modes,
)
from bijux_phylogenetics.comparative.ou_trait_evolution import (
    OUTraitEvolutionSummaryReport,
    summarize_ou_trait_evolution,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.simulation import (
    ContinuousTraitSimulationReport,
    simulate_brownian_traits,
    simulate_early_burst_traits,
    simulate_ou_traits,
    write_continuous_trait_table,
)

_TRAIT_NAME = "value"


@dataclass(slots=True)
class ContinuousModeRecoveryScenario:
    """One deterministic simulation-recovery case for continuous trait models."""

    case_id: str
    label: str
    generating_model: str
    expected_selected_model: str
    root_state: float
    sigma: float
    seed: int
    alpha: float | None = None
    theta: float | None = None
    rate_change: float | None = None
    ou_bounds: tuple[float, float] = (0.0, 10.0)
    early_burst_bounds: tuple[float, float] = (0.0, 50.0)
    parameter_tolerances: dict[str, float] = field(default_factory=dict)
    expected_warning_kinds: list[str] = field(default_factory=list)
    notes: str = ""


@dataclass(slots=True)
class ContinuousModeRecoveryParameterRow:
    """One recovered model parameter compared directly against simulation truth."""

    case_id: str
    generating_model: str
    fitted_model: str
    parameter: str
    true_value: float
    estimated_value: float
    absolute_error: float
    relative_error: float
    tolerance: float
    within_tolerance: bool
    interpretation: str


@dataclass(slots=True)
class ContinuousModeRecoveryModelChoiceRow:
    """One candidate-model row from the BM-versus-OU-versus-EB comparison."""

    case_id: str
    generating_model: str
    expected_selected_model: str
    model: str
    parameter_count: int
    log_likelihood: float
    aic: float
    aicc: float
    selected: bool


@dataclass(slots=True)
class ContinuousModeRecoveryWarningRow:
    """One identifiability or recovery warning emitted during model review."""

    case_id: str
    fitted_model: str
    kind: str
    message: str


@dataclass(slots=True)
class ContinuousModeRecoveryCaseReport:
    """Full recovery review for one deterministic simulation case."""

    scenario: ContinuousModeRecoveryScenario
    traits_path: Path
    simulation: ContinuousTraitSimulationReport
    brownian_fit: BrownianTraitEvolutionSummaryReport
    ou_fit: OUTraitEvolutionSummaryReport
    early_burst_fit: EarlyBurstTraitEvolutionSummaryReport
    model_comparison: ContinuousEvolutionaryModeComparisonReport
    parameter_rows: list[ContinuousModeRecoveryParameterRow]
    model_choice_rows: list[ContinuousModeRecoveryModelChoiceRow]
    warning_rows: list[ContinuousModeRecoveryWarningRow]
    selected_model: str
    selection_matches_expectation: bool
    expected_warning_kinds_present: bool


@dataclass(slots=True)
class ContinuousModeRecoveryReport:
    """Recovery review over one shared tree and multiple deterministic cases."""

    tree_path: Path
    case_reports: list[ContinuousModeRecoveryCaseReport]


def run_continuous_mode_recovery(
    tree_path: Path,
    scenarios: list[ContinuousModeRecoveryScenario],
) -> ContinuousModeRecoveryReport:
    """Simulate, refit, and compare BM, OU, and early-burst cases on one tree."""
    with TemporaryDirectory(prefix="continuous-mode-recovery-") as temporary_root:
        temporary_path = Path(temporary_root)
        case_reports = [
            _run_case(tree_path, scenario, temporary_path)
            for scenario in scenarios
        ]
    return ContinuousModeRecoveryReport(
        tree_path=tree_path,
        case_reports=case_reports,
    )


def write_continuous_mode_recovery_summary_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write one recovery summary row per deterministic simulation case."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "label",
            "generating_model",
            "expected_selected_model",
            "selected_model",
            "selection_matches_expectation",
            "parameter_row_count",
            "parameter_pass_count",
            "expected_warning_count",
            "expected_warning_kinds_present",
            "warning_count",
            "notes",
        ],
        rows=[
            {
                "case_id": case.scenario.case_id,
                "label": case.scenario.label,
                "generating_model": case.scenario.generating_model,
                "expected_selected_model": case.scenario.expected_selected_model,
                "selected_model": case.selected_model,
                "selection_matches_expectation": str(
                    case.selection_matches_expectation
                ).lower(),
                "parameter_row_count": str(len(case.parameter_rows)),
                "parameter_pass_count": str(
                    sum(1 for row in case.parameter_rows if row.within_tolerance)
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
    """Write all parameter-recovery rows across the deterministic cases."""
    rows = [
        {
            "case_id": row.case_id,
            "generating_model": row.generating_model,
            "fitted_model": row.fitted_model,
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
            "fitted_model",
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


def write_continuous_mode_recovery_model_choice_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write the BM/OU/EB model-choice rows for each recovery case."""
    rows = [
        {
            "case_id": row.case_id,
            "generating_model": row.generating_model,
            "expected_selected_model": row.expected_selected_model,
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


def write_continuous_mode_recovery_warning_table(
    path: Path,
    report: ContinuousModeRecoveryReport,
) -> Path:
    """Write the identifiability warnings observed across recovery cases."""
    return write_taxon_rows(
        path,
        columns=["case_id", "fitted_model", "kind", "message"],
        rows=[
            {
                "case_id": row.case_id,
                "fitted_model": row.fitted_model,
                "kind": row.kind,
                "message": row.message,
            }
            for case in report.case_reports
            for row in case.warning_rows
        ],
    )


def _run_case(
    tree_path: Path,
    scenario: ContinuousModeRecoveryScenario,
    temporary_root: Path,
) -> ContinuousModeRecoveryCaseReport:
    simulation = _simulate_case(tree_path, scenario)
    traits_path = write_continuous_trait_table(
        temporary_root / f"{scenario.case_id}-traits.tsv",
        simulation,
    )
    brownian_fit = summarize_brownian_trait_evolution(
        tree_path,
        traits_path,
        trait=_TRAIT_NAME,
    )
    ou_fit = summarize_ou_trait_evolution(
        tree_path,
        traits_path,
        trait=_TRAIT_NAME,
    )
    early_burst_fit = summarize_early_burst_trait_evolution(
        tree_path,
        traits_path,
        trait=_TRAIT_NAME,
        rate_change_bounds=scenario.early_burst_bounds,
    )
    model_comparison = compare_continuous_evolutionary_modes(
        tree_path,
        traits_path,
        trait=_TRAIT_NAME,
        ou_bounds=scenario.ou_bounds,
        early_burst_bounds=scenario.early_burst_bounds,
    )
    parameter_rows = _build_parameter_rows(
        scenario=scenario,
        brownian_fit=brownian_fit,
        ou_fit=ou_fit,
        early_burst_fit=early_burst_fit,
    )
    warning_rows = _build_warning_rows(
        case_id=scenario.case_id,
        ou_fit=ou_fit,
        early_burst_fit=early_burst_fit,
    )
    observed_warning_kinds = {row.kind for row in warning_rows}
    return ContinuousModeRecoveryCaseReport(
        scenario=scenario,
        traits_path=traits_path,
        simulation=simulation,
        brownian_fit=brownian_fit,
        ou_fit=ou_fit,
        early_burst_fit=early_burst_fit,
        model_comparison=model_comparison,
        parameter_rows=parameter_rows,
        model_choice_rows=[
            ContinuousModeRecoveryModelChoiceRow(
                case_id=scenario.case_id,
                generating_model=scenario.generating_model,
                expected_selected_model=scenario.expected_selected_model,
                model=row.model,
                parameter_count=row.parameter_count,
                log_likelihood=row.log_likelihood,
                aic=row.aic,
                aicc=row.aicc,
                selected=row.selected,
            )
            for row in model_comparison.rows
        ],
        warning_rows=warning_rows,
        selected_model=model_comparison.better_model,
        selection_matches_expectation=(
            model_comparison.better_model == scenario.expected_selected_model
        ),
        expected_warning_kinds_present=all(
            kind in observed_warning_kinds for kind in scenario.expected_warning_kinds
        ),
    )


def _simulate_case(
    tree_path: Path,
    scenario: ContinuousModeRecoveryScenario,
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
    raise ValueError(f"unsupported generating_model: {scenario.generating_model}")


def _build_parameter_rows(
    *,
    scenario: ContinuousModeRecoveryScenario,
    brownian_fit: BrownianTraitEvolutionSummaryReport,
    ou_fit: OUTraitEvolutionSummaryReport,
    early_burst_fit: EarlyBurstTraitEvolutionSummaryReport,
) -> list[ContinuousModeRecoveryParameterRow]:
    if not scenario.parameter_tolerances:
        return []
    rows: list[ContinuousModeRecoveryParameterRow] = []
    if scenario.generating_model == "brownian":
        rows.append(
            _parameter_row(
                case_id=scenario.case_id,
                generating_model=scenario.generating_model,
                fitted_model="brownian",
                parameter="sigma_squared",
                true_value=scenario.sigma**2,
                estimated_value=brownian_fit.sigma_squared,
                tolerance=scenario.parameter_tolerances["sigma_squared"],
                interpretation="Brownian sigma squared is recovered by refitting the Brownian model to the simulated tip values.",
            )
        )
        return rows
    if scenario.generating_model == "ornstein-uhlenbeck":
        rows.extend(
            [
                _parameter_row(
                    case_id=scenario.case_id,
                    generating_model=scenario.generating_model,
                    fitted_model="ornstein-uhlenbeck",
                    parameter="alpha",
                    true_value=scenario.alpha or 0.0,
                    estimated_value=ou_fit.alpha,
                    tolerance=scenario.parameter_tolerances["alpha"],
                    interpretation="OU alpha is recovered by refitting the OU model to the simulated tip values.",
                ),
                _parameter_row(
                    case_id=scenario.case_id,
                    generating_model=scenario.generating_model,
                    fitted_model="ornstein-uhlenbeck",
                    parameter="sigma_squared",
                    true_value=scenario.sigma**2,
                    estimated_value=ou_fit.sigma_squared,
                    tolerance=scenario.parameter_tolerances["sigma_squared"],
                    interpretation="OU sigma squared is recovered by refitting the OU model to the simulated tip values.",
                ),
                _parameter_row(
                    case_id=scenario.case_id,
                    generating_model=scenario.generating_model,
                    fitted_model="ornstein-uhlenbeck",
                    parameter="theta",
                    true_value=scenario.theta or 0.0,
                    estimated_value=ou_fit.theta,
                    tolerance=scenario.parameter_tolerances["theta"],
                    interpretation="OU optimum theta is recovered by refitting the OU model to the simulated tip values.",
                ),
            ]
        )
        return rows
    rows.append(
        _parameter_row(
            case_id=scenario.case_id,
            generating_model=scenario.generating_model,
            fitted_model="early-burst",
            parameter="rate_change",
            true_value=scenario.rate_change or 0.0,
            estimated_value=early_burst_fit.rate_change,
            tolerance=scenario.parameter_tolerances["rate_change"],
            interpretation="The early-burst rate-change parameter is recovered by refitting the early-burst model to the simulated tip values.",
        )
    )
    return rows


def _parameter_row(
    *,
    case_id: str,
    generating_model: str,
    fitted_model: str,
    parameter: str,
    true_value: float,
    estimated_value: float,
    tolerance: float,
    interpretation: str,
) -> ContinuousModeRecoveryParameterRow:
    absolute_error = abs(estimated_value - true_value)
    relative_error = 0.0 if true_value == 0.0 else absolute_error / abs(true_value)
    return ContinuousModeRecoveryParameterRow(
        case_id=case_id,
        generating_model=generating_model,
        fitted_model=fitted_model,
        parameter=parameter,
        true_value=true_value,
        estimated_value=estimated_value,
        absolute_error=absolute_error,
        relative_error=relative_error,
        tolerance=tolerance,
        within_tolerance=absolute_error <= tolerance,
        interpretation=interpretation,
    )


def _build_warning_rows(
    *,
    case_id: str,
    ou_fit: OUTraitEvolutionSummaryReport,
    early_burst_fit: EarlyBurstTraitEvolutionSummaryReport,
) -> list[ContinuousModeRecoveryWarningRow]:
    return [
        *[
            ContinuousModeRecoveryWarningRow(
                case_id=case_id,
                fitted_model="ornstein-uhlenbeck",
                kind=warning.kind,
                message=warning.message,
            )
            for warning in ou_fit.identifiability_warnings
        ],
        *[
            ContinuousModeRecoveryWarningRow(
                case_id=case_id,
                fitted_model="early-burst",
                kind=warning.kind,
                message=warning.message,
            )
            for warning in early_burst_fit.identifiability_warnings
        ],
    ]


def _format_number(value: float) -> str:
    return format(value, ".15g")
