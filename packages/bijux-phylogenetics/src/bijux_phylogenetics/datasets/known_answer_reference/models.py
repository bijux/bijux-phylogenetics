from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.continuous import ContinuousAncestralReport
from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport
from bijux_phylogenetics.comparative.continuous import (
    BrownianTraitEvolutionSummaryReport,
    OUTraitEvolutionSummaryReport,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    DiscreteStateEvolutionReport,
)
from bijux_phylogenetics.compare.topology import TreeComparisonReport
from bijux_phylogenetics.distance import DistanceTreeBuildReport
from bijux_phylogenetics.ecology import HostSwitchingReport

DATASET_ID = "known_answer_reference_panel"
DATASET_LABEL = "Known-answer simulation reference panel"
SEQUENCE_TYPE = "dna"
DISTANCE_METHOD = "neighbor-joining"
DISTANCE_MODEL = "p-distance"
CONTINUOUS_TRAIT = "value"
DISCRETE_TRAIT = "state"


@dataclass(slots=True)
class KnownAnswerContinuousNodeTruth:
    """One true continuous node value stored with the packaged simulation panel."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    true_value: float


@dataclass(slots=True)
class KnownAnswerDiscreteNodeTruth:
    """One true discrete node state stored with the packaged simulation panel."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    true_state: str


@dataclass(slots=True)
class KnownAnswerTransitionTruth:
    """One true simulated branch-change row stored with the packaged panel."""

    parent_node: str
    child_node: str
    branch_length: float
    source_state: str
    target_state: str
    changed: bool
    event_count: int


@dataclass(slots=True)
class KnownAnswerRecoveryThreshold:
    """One declared pass or fail threshold for the known-answer suite."""

    metric: str
    comparator: str
    threshold: str
    rationale: str


@dataclass(slots=True)
class KnownAnswerParameterRecoveryRow:
    """One recovered continuous parameter compared directly against truth."""

    parameter: str
    true_value: float
    estimated_value: float
    absolute_error: float
    relative_error: float
    interpretation: str


@dataclass(slots=True)
class KnownAnswerContinuousNodeRecoveryRow:
    """One internal-node continuous ancestral estimate compared against truth."""

    node: str
    descendant_taxa: list[str]
    true_value: float
    estimated_value: float
    absolute_error: float
    standard_error: float
    lower_95_interval: float
    upper_95_interval: float
    confidence: float


@dataclass(slots=True)
class KnownAnswerDiscreteNodeRecoveryRow:
    """One internal-node discrete ancestral estimate compared against truth."""

    node: str
    descendant_taxa: list[str]
    true_state: str
    estimated_state: str
    true_state_probability: float
    confidence: float
    correct: bool
    ambiguous: bool


@dataclass(slots=True)
class KnownAnswerTransitionRecoveryRow:
    """One branchwise transition recovery row compared against stored truth."""

    parent_node: str
    child_node: str
    true_transition: str
    estimated_transition: str
    true_changed: bool
    estimated_changed: bool
    true_event_count: int
    estimated_event_count: int
    correct: bool


@dataclass(slots=True)
class KnownAnswerThresholdEvaluationRow:
    """One evaluated known-answer threshold against observed recovery metrics."""

    metric: str
    comparator: str
    threshold: str
    observed_value: str
    passed: bool
    rationale: str


@dataclass(slots=True)
class KnownAnswerReferenceDataset:
    """Packaged deterministic simulation panel with stored truth artifacts."""

    dataset_id: str
    label: str
    dataset_root: Path
    true_tree_path: Path
    alignment_path: Path
    continuous_traits_path: Path
    ou_traits_path: Path
    discrete_traits_path: Path
    host_traits_path: Path
    geographic_traits_path: Path
    true_parameters_path: Path
    true_continuous_nodes_path: Path
    true_ou_nodes_path: Path
    true_discrete_nodes_path: Path
    true_host_nodes_path: Path
    true_geographic_nodes_path: Path
    true_host_switch_events_path: Path
    true_geographic_transition_events_path: Path
    recovery_thresholds_path: Path
    reference_output_root: Path
    taxon_count: int
    sequence_length: int
    sequence_type: str
    distance_method: str
    distance_model: str
    source_summary: str


@dataclass(slots=True)
class KnownAnswerReferenceExportResult:
    """Materialized copy of the packaged simulation dataset."""

    output_root: Path
    readme_path: Path
    true_tree_path: Path
    alignment_path: Path
    continuous_traits_path: Path
    ou_traits_path: Path
    discrete_traits_path: Path
    host_traits_path: Path
    geographic_traits_path: Path
    true_parameters_path: Path
    true_continuous_nodes_path: Path
    true_ou_nodes_path: Path
    true_discrete_nodes_path: Path
    true_host_nodes_path: Path
    true_geographic_nodes_path: Path
    true_host_switch_events_path: Path
    true_geographic_transition_events_path: Path
    recovery_thresholds_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class KnownAnswerReferenceWorkflowReport:
    """Recovery workflow run over the packaged known-answer simulation panel."""

    dataset: KnownAnswerReferenceDataset
    distance_tree_build: DistanceTreeBuildReport
    distance_tree_newick: str
    tree_recovery: TreeComparisonReport
    brownian_fit: BrownianTraitEvolutionSummaryReport
    ou_fit: OUTraitEvolutionSummaryReport
    continuous_ancestral: ContinuousAncestralReport
    discrete_ancestral: DiscreteAncestralReport
    host_switching: HostSwitchingReport
    geographic_states: DiscreteStateEvolutionReport
    parameter_recovery_rows: list[KnownAnswerParameterRecoveryRow]
    continuous_node_recovery_rows: list[KnownAnswerContinuousNodeRecoveryRow]
    discrete_node_recovery_rows: list[KnownAnswerDiscreteNodeRecoveryRow]
    host_node_recovery_rows: list[KnownAnswerDiscreteNodeRecoveryRow]
    host_event_recovery_rows: list[KnownAnswerTransitionRecoveryRow]
    geographic_node_recovery_rows: list[KnownAnswerDiscreteNodeRecoveryRow]
    geographic_event_recovery_rows: list[KnownAnswerTransitionRecoveryRow]
    threshold_evaluation_rows: list[KnownAnswerThresholdEvaluationRow]


@dataclass(slots=True)
class KnownAnswerReferenceWorkflowBundle:
    """Written recovery outputs for the packaged known-answer simulation panel."""

    output_root: Path
    rooted_topology_equal: bool
    same_unrooted_topology: bool
    same_taxa_different_rooting: bool
    robinson_foulds_distance: int
    continuous_internal_node_mean_absolute_error: float
    discrete_internal_node_accuracy: float
    discrete_mean_true_state_probability: float
    host_internal_node_accuracy: float
    host_event_accuracy: float
    geographic_internal_node_accuracy: float
    geographic_event_accuracy: float
    parameter_row_count: int
    threshold_pass_count: int
    threshold_row_count: int
    workflow_summary_path: Path
    distance_tree_path: Path
    tree_recovery_path: Path
    parameter_recovery_path: Path
    brownian_fit_summary_path: Path
    ou_fit_summary_path: Path
    continuous_ancestral_summary_path: Path
    continuous_ancestral_uncertainty_path: Path
    continuous_node_recovery_path: Path
    discrete_ancestral_summary_path: Path
    discrete_ancestral_probability_path: Path
    discrete_node_recovery_path: Path
    host_switch_summary_path: Path
    host_state_nodes_path: Path
    host_switch_branches_path: Path
    host_node_recovery_path: Path
    host_event_recovery_path: Path
    geographic_ancestral_summary_path: Path
    geographic_state_probability_path: Path
    geographic_transition_summary_path: Path
    geographic_node_recovery_path: Path
    geographic_event_recovery_path: Path
    threshold_evaluation_path: Path


@dataclass(slots=True)
class KnownAnswerReferenceDemoResult:
    """Dataset export plus recovery workflow outputs for the public simulation demo."""

    output_root: Path
    dataset: KnownAnswerReferenceDataset
    dataset_export: KnownAnswerReferenceExportResult
    workflow_bundle: KnownAnswerReferenceWorkflowBundle
    overview_path: Path
