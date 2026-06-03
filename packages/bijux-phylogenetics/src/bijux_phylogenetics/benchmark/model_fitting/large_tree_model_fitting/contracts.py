from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class LargeTreeModelFittingThreshold:
    """Performance budget for one governed large-tree model-fitting benchmark case."""

    max_runtime_seconds: float
    max_peak_memory_bytes: int
    max_optimizer_step_count: int


@dataclass(slots=True)
class LargeTreeModelFittingObservation:
    """One owned large-tree model-fitting benchmark result with geiger comparison."""

    case_id: str
    tier: str
    trait_kind: str
    fit_surface: str
    taxon_count: int
    status: str
    runtime_seconds: float | None
    peak_memory_bytes: int | None
    memory_observation_kind: str | None
    optimizer_name: str | None
    optimizer_iteration_count: int | None
    optimizer_function_evaluation_count: int | None
    converged: bool | None
    hit_lower_parameter_boundary: bool | None
    hit_upper_parameter_boundary: bool | None
    unstable_review: bool
    too_slow_review: bool
    stable_conclusion_supported: bool | None
    threshold: LargeTreeModelFittingThreshold
    runtime_within_threshold: bool | None
    peak_memory_within_threshold: bool | None
    optimizer_step_within_threshold: bool | None
    performance_threshold_passed: bool | None
    geiger_reference_available: bool
    geiger_runtime_seconds: float | None
    geiger_optimizer_step_count: int | None
    geiger_parameter_name: str | None
    geiger_parameter_value: float | None
    geiger_rate: float | None
    geiger_log_likelihood: float | None
    geiger_aic: float | None
    geiger_aicc: float | None
    parameter_delta: float | None
    rate_delta: float | None
    log_likelihood_delta: float | None
    aic_delta: float | None
    geiger_match_tolerance: float | None
    matches_geiger_reference: bool | None
    notes: list[str]


@dataclass(slots=True)
class LargeTreeModelFittingBenchmarkReport:
    """Tiered large-tree model-fitting benchmark with governed geiger comparison."""

    tier: str
    observations: list[LargeTreeModelFittingObservation]
    case_count: int
    geiger_match_case_count: int
    threshold_pass_case_count: int
    too_slow_case_count: int
    unstable_case_count: int
    limitations: list[str]


@dataclass(slots=True)
class LargeTreeModelFittingBenchmarkBundle:
    """Written artifact bundle for one governed large-tree model-fitting benchmark tier."""

    output_root: Path
    summary_path: Path
    observation_table_path: Path
