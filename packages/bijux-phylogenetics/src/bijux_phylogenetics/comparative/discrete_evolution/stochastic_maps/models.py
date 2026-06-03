from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class StochasticMapTransitionEvent:
    branch_index: int
    parent_node: str
    child_node: str
    source_state: str
    target_state: str
    event_time_fraction: float


@dataclass(slots=True)
class StochasticMapStateSegment:
    branch_index: int
    parent_node: str
    child_node: str
    state: str
    start_time_fraction: float
    end_time_fraction: float
    duration: float


@dataclass(slots=True)
class StochasticMapBranchHistory:
    branch_index: int
    parent_node: str
    child_node: str
    branch_length: float
    start_state: str
    end_state: str
    event_count: int
    events: list[StochasticMapTransitionEvent]
    segments: list[StochasticMapStateSegment]


@dataclass(slots=True)
class StochasticMapReplicate:
    replicate_index: int
    root_state: str
    total_transition_count: int
    transition_counts: dict[str, int]
    state_time_totals: dict[str, float]
    branch_histories: list[StochasticMapBranchHistory]


@dataclass(slots=True)
class StochasticMapSummaryRow:
    transition: str
    mean_count: float
    lower_95_interval: float
    upper_95_interval: float
    minimum_count: int
    maximum_count: int
    presence_fraction: float


@dataclass(slots=True)
class StochasticMapStateTimeRow:
    state: str
    mean_time: float
    lower_95_interval: float
    upper_95_interval: float
    minimum_time: float
    maximum_time: float


@dataclass(slots=True)
class StochasticMapBranchOccupancyRow:
    branch_index: int
    parent_node: str
    child_node: str
    state: str
    branch_length: float
    mean_time: float
    lower_95_interval: float
    upper_95_interval: float
    minimum_time: float
    maximum_time: float
    mean_fraction: float
    presence_fraction: float


@dataclass(slots=True)
class StochasticMapTransitionCountMatrixRow:
    replicate_index: int
    total_transition_count: int
    transition_counts: dict[str, int]


@dataclass(slots=True)
class StochasticMapBranchTransitionCountRow:
    branch_index: int
    parent_node: str
    child_node: str
    transition: str
    mean_count: float
    lower_95_interval: float
    upper_95_interval: float
    minimum_count: int
    maximum_count: int
    presence_fraction: float


@dataclass(slots=True)
class StochasticMapSimulationFailure:
    replicate_index: int
    branch_index: int
    parent_node: str
    child_node: str
    source_state: str
    target_state: str
    branch_length: float
    attempt_count: int
    reason: str


@dataclass(slots=True)
class StochasticMapModelFitAudit:
    state_order: list[str]
    allowed_transitions: list[str]
    parameter_count: int
    log_likelihood: float
    aic: float
    aicc: float
    overparameterized: bool
    optimizer_converged: bool
    optimizer_iteration_count: int
    optimizer_function_evaluation_count: int
    optimizer_hit_lower_parameter_bound: bool
    optimizer_hit_upper_parameter_bound: bool
    baseline_model: str | None
    baseline_aic: float | None
    baseline_delta_aic: float | None
    preferred_model_by_aic: str | None
    warnings: list[str]


@dataclass(slots=True)
class StochasticMapSummaryReport:
    replicate_count: int
    mean_total_transition_count: float
    lower_95_total_transition_count: float
    upper_95_total_transition_count: float
    rows: list[StochasticMapSummaryRow]
    state_time_rows: list[StochasticMapStateTimeRow]
    branch_occupancy_rows: list[StochasticMapBranchOccupancyRow]
    simulation_failure_count: int
    warnings: list[str]


@dataclass(slots=True)
class StochasticMapCollectionReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    state_ordering: str
    ordered_states: list[str]
    replicates: int
    seed: int
    conditioned_on_node_estimates: bool
    fit_audit: StochasticMapModelFitAudit
    warnings: list[str]
    maps: list[StochasticMapReplicate]
    failures: list[StochasticMapSimulationFailure]
    summary: StochasticMapSummaryReport


@dataclass(slots=True)
class StochasticMapTransitionCountReport:
    replicate_count: int
    mean_total_transition_count: float
    lower_95_total_transition_count: float
    upper_95_total_transition_count: float
    transition_order: list[str]
    matrix_rows: list[StochasticMapTransitionCountMatrixRow]
    aggregate_rows: list[StochasticMapSummaryRow]
    branch_rows: list[StochasticMapBranchTransitionCountRow]
    warnings: list[str]


@dataclass(slots=True)
class StochasticMapBranchProbabilityRow:
    branch_index: int
    parent_node: str
    child_node: str
    state: str
    branch_length: float
    mean_probability: float
    lower_95_probability: float
    upper_95_probability: float
    minimum_probability: float
    maximum_probability: float
    presence_fraction: float


@dataclass(slots=True)
class StochasticMapDensitySliceRow:
    branch_index: int
    parent_node: str
    child_node: str
    branch_length: float
    slice_index: int
    start_depth: float
    end_depth: float
    start_time_fraction: float
    end_time_fraction: float
    posterior_probability: float
    posterior_uncertainty: float


@dataclass(slots=True)
class StochasticMapDensityBranchRow:
    branch_index: int
    parent_node: str
    child_node: str
    branch_length: float
    focal_state: str
    baseline_state: str | None
    mean_posterior_probability: float
    minimum_posterior_probability: float
    maximum_posterior_probability: float
    uncertainty: float
    slice_count: int


@dataclass(slots=True)
class StochasticMapDensityReport:
    replicate_count: int
    resolution: int
    total_tree_depth: float
    state_order: list[str]
    focal_state: str | None
    baseline_state: str | None
    branch_state_rows: list[StochasticMapBranchProbabilityRow]
    density_rows: list[StochasticMapDensitySliceRow]
    branch_rows: list[StochasticMapDensityBranchRow]
    warnings: list[str]


@dataclass(slots=True)
class StochasticMapDensityArtifactResult:
    output_path: Path
    svg_path: Path
    format: str
    layout: str
    focal_state: str
    baseline_state: str | None
    branch_count: int
    rendered_branch_color_count: int
