from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from bijux_phylogenetics.simulation import (
    DiscreteHistoryRateRow,
    DiscreteTraitSimulationReport,
)

RECOVERY_TRAIT_NAME = "state"
DEFAULT_CANDIDATE_MODELS = (
    "equal-rates",
    "symmetric",
    "all-rates-different",
)


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
    transform: str | None = None
    transform_parameter_value: float | None = None
    candidate_models: tuple[str, ...] = DEFAULT_CANDIDATE_MODELS
    rate_tolerance: float | None = None
    parameter_tolerances: dict[str, float] = field(default_factory=dict)
    lambda_bounds: tuple[float, float] = (0.0, 1.0)
    kappa_bounds: tuple[float, float] = (0.0, 1.0)
    delta_bounds: tuple[float, float] = (0.006737947, 3.0)
    early_burst_bounds: tuple[float, float] = (-10.0, 10.0)
    expected_overparameterized: bool = False
    expected_warning_kinds: list[str] = field(default_factory=list)
    notes: str = ""


@dataclass(slots=True)
class DiscreteModeRecoveryParameterRow:
    """One truth-versus-fit transform-parameter recovery row for one engine."""

    case_id: str
    generating_model: str
    transform: str
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
class DiscreteModeRecoveryParameterComparisonRow:
    """One paired Bijux-versus-geiger transform-parameter recovery comparison."""

    case_id: str
    generating_model: str
    transform: str
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
class DiscreteModeRecoveryRateRow:
    """One truth-versus-fit transition-rate recovery row for one engine."""

    case_id: str
    generating_model: str
    transform: str
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
    transform: str
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
    transform: str
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
    parameter_rows: list[DiscreteModeRecoveryParameterRow]
    parameter_comparison_rows: list[DiscreteModeRecoveryParameterComparisonRow]
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
    parameter_name: str | None
    parameter_value: float | None
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
