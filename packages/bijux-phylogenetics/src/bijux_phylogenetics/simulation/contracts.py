from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from bijux_phylogenetics.phylo.alignment import AlignmentRecord


@dataclass(frozen=True, slots=True)
class SimulatedTreeRecord:
    index: int
    newick: str
    tree_height_branch_length: float
    total_branch_length: float
    mean_branch_length: float
    median_branch_length: float
    minimum_branch_length: float
    maximum_branch_length: float
    cherry_count: int
    sackin_imbalance_index: int
    normalized_colless_imbalance: float


@dataclass(frozen=True, slots=True)
class TreeSimulationEnvelopeMetric:
    metric: str
    sample_scope: str
    observation_count: int
    mean: float
    standard_deviation: float
    minimum: float
    median: float
    maximum: float


@dataclass(frozen=True, slots=True)
class CoalescentWaitingTimeSummaryRow:
    lineage_count: int
    coalescent_rate: float
    expected_waiting_time: float
    observation_count: int
    mean_waiting_time: float
    standard_deviation: float
    minimum_waiting_time: float
    median_waiting_time: float
    maximum_waiting_time: float
    absolute_error: float
    relative_error: float
    within_tolerance: bool


@dataclass(frozen=True, slots=True)
class CoalescentSkylineSummaryRow:
    interval: str
    lineage_count: int
    duration: float
    effective_population_size_estimate: float
    observation_count: int
    relative_error: float
    uncertainty_flag: str


@dataclass(slots=True)
class TreeSimulationReport:
    model: str
    tree_count: int
    tip_count: int
    seed: int
    records: list[SimulatedTreeRecord]
    branch_length_model: str | None = None
    birth_rate: float | None = None
    death_rate: float | None = None
    population_size: float | None = None
    rooted: bool = True
    binary: bool = True
    pooled_branch_count: int = 0
    envelope_metrics: list[TreeSimulationEnvelopeMetric] = field(default_factory=list)
    coalescent_waiting_time_tolerance: float | None = None
    coalescent_waiting_time_rows: list[CoalescentWaitingTimeSummaryRow] = field(
        default_factory=list
    )
    coalescent_skyline_rows: list[CoalescentSkylineSummaryRow] = field(
        default_factory=list
    )


@dataclass(frozen=True, slots=True)
class MultispeciesCoalescentSampleRow:
    species_taxon: str
    sample_count: int
    gene_taxa: list[str]


@dataclass(frozen=True, slots=True)
class MultispeciesCoalescentEventRow:
    event_index: int
    species_branch: str
    branch_role: str
    descendant_species: list[str]
    population_size: float
    branch_start_age: float
    branch_end_age: float | None
    event_age: float
    waiting_time: float
    input_lineage_count: int
    output_lineage_count: int
    left_gene_clade: list[str]
    right_gene_clade: list[str]
    resulting_gene_clade: list[str]


@dataclass(frozen=True, slots=True)
class MultispeciesCoalescentBranchRow:
    species_branch: str
    branch_role: str
    descendant_species: list[str]
    branch_duration: float | None
    population_size: float
    lineage_count_entering: int
    coalescent_event_count: int
    lineage_count_exiting: int
    extra_lineage_count: int
    included_in_deep_coalescence_total: bool


@dataclass(slots=True)
class MultispeciesCoalescentReport:
    model: str
    species_tree_path: Path
    seed: int
    species_tip_count: int
    gene_tip_count: int
    default_population_size: float
    deep_coalescence_total: int
    sample_rows: list[MultispeciesCoalescentSampleRow]
    branch_rows: list[MultispeciesCoalescentBranchRow]
    event_rows: list[MultispeciesCoalescentEventRow]


@dataclass(frozen=True, slots=True)
class SimulatedContinuousTrait:
    taxon: str
    value: float


@dataclass(frozen=True, slots=True)
class SimulatedContinuousNode:
    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    value: float


@dataclass(slots=True)
class ContinuousTraitSimulationReport:
    model: str
    tree_path: Path
    tip_count: int
    seed: int
    root_state: float
    sigma: float
    sigma_squared: float
    alpha: float | None
    theta: float | None
    rate_change: float | None
    traits: list[SimulatedContinuousTrait]
    node_values: list[SimulatedContinuousNode]


@dataclass(frozen=True, slots=True)
class ContinuousTraitSimulationSummaryRow:
    row_kind: str
    label: str
    mean_value: float | None = None
    standard_deviation: float | None = None
    minimum: float | None = None
    median: float | None = None
    maximum: float | None = None
    covariance: float | None = None
    correlation: float | None = None


@dataclass(slots=True)
class ContinuousTraitSimulationCollectionReport:
    model: str
    tree_path: Path
    tip_count: int
    branch_count: int
    replicate_count: int
    seed: int
    root_state: float
    sigma: float
    sigma_squared: float
    simulations: list[ContinuousTraitSimulationReport]
    rows: list[ContinuousTraitSimulationSummaryRow]


@dataclass(frozen=True, slots=True)
class SimulatedCorrelatedContinuousTrait:
    taxon: str
    trait: str
    value: float


@dataclass(slots=True)
class CorrelatedContinuousTraitSimulationReport:
    model: str
    tree_path: Path
    tip_count: int
    trait_names: list[str]
    seed: int
    root_states: list[float]
    evolutionary_covariance_matrix: list[list[float]]
    traits: list[SimulatedCorrelatedContinuousTrait]


@dataclass(slots=True)
class CorrelatedContinuousTraitSimulationCollectionReport:
    model: str
    tree_path: Path
    tip_count: int
    branch_count: int
    trait_names: list[str]
    replicate_count: int
    seed: int
    root_states: list[float]
    evolutionary_covariance_matrix: list[list[float]]
    simulations: list[CorrelatedContinuousTraitSimulationReport]
    rows: list[ContinuousTraitSimulationSummaryRow]


@dataclass(frozen=True, slots=True)
class SimulatedDiscreteTrait:
    taxon: str
    state: str


@dataclass(frozen=True, slots=True)
class SimulatedDiscreteNode:
    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    state: str


@dataclass(frozen=True, slots=True)
class SimulatedDiscreteTransitionEvent:
    parent_node: str
    child_node: str
    source_state: str
    target_state: str
    event_index: int
    branch_distance: float = 0.0


@dataclass(frozen=True, slots=True)
class SimulatedDiscreteStateSegment:
    parent_node: str
    child_node: str
    state: str
    start_distance: float
    end_distance: float
    duration: float


@dataclass(frozen=True, slots=True)
class DiscreteHistoryRateRow:
    source_state: str
    target_state: str
    rate: float


@dataclass(frozen=True, slots=True)
class SimulatedDiscreteBranchHistory:
    parent_node: str
    child_node: str
    branch_length: float
    start_state: str
    end_state: str
    changed: bool
    event_count: int
    events: list[SimulatedDiscreteTransitionEvent]
    segments: list[SimulatedDiscreteStateSegment]


@dataclass(slots=True)
class DiscreteTraitSimulationReport:
    model: str
    tree_path: Path
    tip_count: int
    seed: int
    states: list[str]
    transition_rate: float | None
    root_state: str
    root_state_probabilities: dict[str, float]
    traits: list[SimulatedDiscreteTrait]
    node_states: list[SimulatedDiscreteNode]
    branch_histories: list[SimulatedDiscreteBranchHistory]
    rate_rows: list[DiscreteHistoryRateRow] = field(default_factory=list)
    transform_name: str | None = None
    transform_parameter_name: str | None = None
    transform_parameter_value: float | None = None


@dataclass(frozen=True, slots=True)
class DiscreteHistorySummaryRow:
    row_kind: str
    label: str
    mean_value: float
    lower_95_interval: float
    upper_95_interval: float
    presence_fraction: float


@dataclass(slots=True)
class DiscreteHistorySimulationCollectionReport:
    model: str
    tree_path: Path
    tip_count: int
    branch_count: int
    replicate_count: int
    seed: int
    states: list[str]
    fixed_root_state: str | None
    root_state_probabilities: dict[str, float]
    rate_rows: list[DiscreteHistoryRateRow]
    simulations: list[DiscreteTraitSimulationReport]
    mean_total_transition_count: float
    lower_95_total_transition_count: float
    upper_95_total_transition_count: float
    rows: list[DiscreteHistorySummaryRow]
    transform_name: str | None = None
    transform_parameter_name: str | None = None
    transform_parameter_value: float | None = None


@dataclass(slots=True)
class AlignmentSimulationReport:
    model: str
    tree_path: Path
    tip_count: int
    seed: int
    sequence_length: int
    substitution_rate: float
    inferred_alphabet: str
    records: list[AlignmentRecord]
