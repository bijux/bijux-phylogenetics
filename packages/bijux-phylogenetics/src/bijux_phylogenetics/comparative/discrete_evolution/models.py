from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class StateCodingIssue:
    taxon: str
    raw_state: str
    code: str
    message: str


@dataclass(slots=True)
class StateCodingValidationReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    allowed_states: list[str]
    state_ordering: str
    ordered_states: list[str]
    valid: bool
    issues: list[StateCodingIssue]
    observed_states: list[str]
    usable_taxa: list[str]


@dataclass(slots=True)
class StateCodingAuditRow:
    taxon: str
    raw_state: str
    normalized_state: str | None
    in_tree: bool
    included: bool
    issue_code: str | None
    note: str


@dataclass(slots=True)
class StateCodingAuditReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    state_ordering: str
    ordered_states: list[str]
    coding_map: dict[str, str]
    row_count: int
    included_row_count: int
    excluded_row_count: int
    rows: list[StateCodingAuditRow]


@dataclass(slots=True)
class StateImbalanceWarning:
    code: str
    message: str
    affected_states: list[str]


@dataclass(slots=True)
class StateImbalanceReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    taxon_count: int
    observed_states: list[str]
    state_counts: dict[str, int]
    warnings: list[StateImbalanceWarning]


@dataclass(slots=True)
class TransitionRateRow:
    source_state: str
    target_rates: dict[str, float]


@dataclass(slots=True)
class TransitionRateUncertaintyRow:
    source_state: str
    target_state: str
    estimate: float
    lower_95_interval: float
    upper_95_interval: float
    effective_transition_count: float


@dataclass(slots=True)
class TransitionRateUncertaintyReport:
    model: str
    state_ordering: str
    rows: list[TransitionRateUncertaintyRow]


@dataclass(slots=True)
class SparseStateInstabilityReport:
    sparse_states: list[str]
    zero_support_transitions: list[str]
    warning_count: int
    unstable: bool


@dataclass(slots=True)
class DominantStateBiasReport:
    dominant_states: list[str]
    dominant_fraction: float
    biased: bool
    message: str | None


@dataclass(slots=True)
class GeographicAnalysisReadinessReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    valid: bool
    blockers: list[str]
    warnings: list[str]
    state_ordering: str
    ordered_states: list[str]
    coding_validation: StateCodingValidationReport
    imbalance: StateImbalanceReport
    dominant_state_bias: DominantStateBiasReport
    tree_validation_decision: str


@dataclass(slots=True)
class TransitionModelReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    likelihood_method: str
    state_ordering: str
    ordered_states: list[str]
    state_order: list[str]
    parameter_count: int
    pseudo_log_likelihood: float
    aic: float
    stationary_frequencies: dict[str, float]
    transition_matrix: list[TransitionRateRow]
    uncertainty: TransitionRateUncertaintyReport
    root_state_probabilities: dict[str, float]


@dataclass(slots=True)
class NodeStateEstimate:
    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    most_likely_state: str
    state_probabilities: dict[str, float]
    ambiguous: bool


@dataclass(slots=True)
class TransitionEvent:
    parent_node: str
    child_node: str
    source_state: str
    target_state: str
    changed: bool


@dataclass(slots=True)
class TransitionSupportRow:
    parent_node: str
    child_node: str
    inferred_transition: str
    support: float
    strongly_supported: bool


@dataclass(slots=True)
class TransitionSummaryReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    branch_count: int
    transition_count: int
    strongly_supported_transition_count: int
    transition_counts: dict[str, int]
    strongly_supported_transition_counts: dict[str, int]
    support_rows: list[TransitionSupportRow]
    events: list[TransitionEvent]


@dataclass(slots=True)
class DiscreteStateEvolutionReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    likelihood_method: str
    state_ordering: str
    ordered_states: list[str]
    analysis_tree_newick: str
    taxon_count: int
    observed_states: list[str]
    state_counts: dict[str, int]
    coding_validation: StateCodingValidationReport
    imbalance: StateImbalanceReport
    instability: SparseStateInstabilityReport
    dominant_state_bias: DominantStateBiasReport
    transition_model: TransitionModelReport
    estimates: list[NodeStateEstimate]
    transition_summary: TransitionSummaryReport
    warnings: list[str]


@dataclass(slots=True)
class DiscreteModelComparisonRow:
    model: str
    parameter_count: int
    pseudo_log_likelihood: float
    aic: float
    transition_count: int


@dataclass(slots=True)
class ModelSensitiveRegionRow:
    node: str
    descendant_taxa: list[str]
    left_state: str
    right_state: str
    sensitivity_score: float


@dataclass(slots=True)
class NodeStateDifference:
    node: str
    descendant_taxa: list[str]
    left_state: str
    right_state: str
    differs: bool
    left_probabilities: dict[str, float]
    right_probabilities: dict[str, float]


@dataclass(slots=True)
class DiscreteModelComparisonReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    left_model: str
    right_model: str
    better_model: str
    rows: list[DiscreteModelComparisonRow]
    node_differences: list[NodeStateDifference]
    sensitive_region_count: int
    sensitive_regions: list[ModelSensitiveRegionRow]


@dataclass(slots=True)
class DiscreteEvolutionNarrative:
    summary: str
    transition_summary: str
    interpretation_boundary: str
    caveats: list[str]


@dataclass(slots=True)
class BiogeographicComputedResult:
    label: str
    value: str


@dataclass(slots=True)
class BiogeographicInterpretationReport:
    tree_path: Path
    traits_path: Path
    trait: str
    model: str
    compare_model: str | None
    computed_results: list[BiogeographicComputedResult]
    model_sensitive_regions: list[ModelSensitiveRegionRow]
    coding_audit_summary: dict[str, int]
    readiness_blockers: list[str]
    caveats: list[str]
    interpretation_guidance: list[str]


@dataclass(slots=True)
class DiscreteEvolutionReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    tree_path: Path
    traits_path: Path
    trait: str
    model: str
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class DiscreteTransitionReferenceRate:
    source_state: str
    target_state: str
    expected_rate: float
    observed_rate: float
    absolute_delta: float


@dataclass(slots=True)
class DiscreteTransitionReferenceObservation:
    label: str
    model: str
    expected_parameter_count: int
    observed_parameter_count: int
    expected_transition_count: int
    observed_transition_count: int
    expected_root_state: str
    observed_root_state: str
    expected_pseudo_log_likelihood: float
    observed_pseudo_log_likelihood: float
    max_rate_delta: float
    rate_rows: list[DiscreteTransitionReferenceRate]
    passed: bool


@dataclass(slots=True)
class DiscreteTransitionReferenceValidationReport:
    case_count: int
    all_passed: bool
    tolerance: float
    observations: list[DiscreteTransitionReferenceObservation]
