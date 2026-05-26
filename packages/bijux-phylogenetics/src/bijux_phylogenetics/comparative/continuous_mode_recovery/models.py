from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from bijux_phylogenetics.simulation import ContinuousTraitSimulationReport

RECOVERY_TRAIT_NAME = "value"
LEGACY_CANDIDATE_MODES = ("brownian", "ornstein-uhlenbeck", "early-burst")
FULL_CANDIDATE_MODES = (
    "brownian",
    "white-noise",
    "pagel-lambda",
    "pagel-kappa",
    "pagel-delta",
    "ornstein-uhlenbeck",
    "early-burst",
)
GEIGER_MODEL_NAMES = {
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
    candidate_modes: tuple[str, ...] = LEGACY_CANDIDATE_MODES
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
    root_state: float | None
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
