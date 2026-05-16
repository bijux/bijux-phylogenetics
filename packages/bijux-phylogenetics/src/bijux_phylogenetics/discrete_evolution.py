from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import random
import tempfile

import numpy

from bijux_phylogenetics.ancestral.common import (
    load_discrete_dataset,
    node_descendant_taxa,
    node_signature,
)
from bijux_phylogenetics.ancestral.discrete import (
    _branch_length,
    _resolve_discrete_model_name,
    _resolve_root_prior,
    _transition_probability_matrix,
)
from bijux_phylogenetics.comparative.discrete_mk import (
    DiscreteMkFitReport,
    fit_discrete_mk_model_from_dataset,
)
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.core.traits import load_tsv_summary
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.errors import AncestralReconstructionError
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.html import write_html_report
from bijux_phylogenetics.render.svg import TreeRenderResult, render_tree_svg

_DEFAULT_ALLOWED_STATE_PATTERN_BLOCKLIST = ("|", "/", ";")
_DEFAULT_STATE_COLORS = (
    "#0f766e",
    "#1d4ed8",
    "#c2410c",
    "#7c3aed",
    "#b91c1c",
    "#047857",
    "#a16207",
    "#0f172a",
)


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


def _normalize_probabilities(scores: dict[str, float]) -> dict[str, float]:
    total = sum(scores.values())
    if total <= 0.0:
        uniform = 1.0 / max(len(scores), 1)
        return {state: float(format(uniform, ".15g")) for state in sorted(scores)}
    return {
        state: float(format(scores[state] / total, ".15g")) for state in sorted(scores)
    }


def _resolve_state_order(
    observed_states: list[str],
    *,
    allowed_states: list[str] | None,
    ordered_states: list[str] | None,
    state_ordering: str,
) -> list[str]:
    if state_ordering == "ordered":
        if ordered_states:
            return ordered_states
        if allowed_states:
            return allowed_states
    return sorted(observed_states)


def _transition_allowed(
    source_state: str,
    target_state: str,
    state_order: list[str],
    *,
    state_ordering: str,
) -> bool:
    if source_state == target_state or state_ordering != "ordered":
        return True
    return abs(state_order.index(source_state) - state_order.index(target_state)) == 1


def _state_support(
    node, states_by_taxon: dict[str, str], state_order: list[str]
) -> dict[str, int]:
    counts = dict.fromkeys(state_order, 0)
    for taxon in node_descendant_taxa(node):
        state = states_by_taxon.get(taxon)
        if state is not None:
            counts[state] += 1
    return counts


def _fitch_candidate_sets(
    tree, states_by_taxon: dict[str, str]
) -> dict[str, list[str]]:
    def downpass(node) -> set[str]:
        if node.is_leaf():
            return {states_by_taxon[node.name]}
        child_sets = [downpass(child) for child in node.children]
        candidate = set(child_sets[0])
        for child_set in child_sets[1:]:
            intersection = candidate & child_set
            if intersection:
                candidate = intersection
            else:
                candidate |= child_set
        return candidate

    return {node_signature(node): sorted(downpass(node)) for node in tree.iter_nodes()}


def _best_supported_state(
    candidate_states: list[str],
    support_counts: dict[str, int],
    priority_weights: dict[str, float],
) -> str:
    return max(
        sorted(candidate_states),
        key=lambda state: (
            support_counts.get(state, 0),
            priority_weights.get(state, 0.0),
            state,
        ),
    )


def _resolve_er_states(
    tree,
    candidate_sets: dict[str, list[str]],
    states_by_taxon: dict[str, str],
    state_order: list[str],
) -> dict[str, str]:
    support_by_node = {
        node_signature(node): _state_support(node, states_by_taxon, state_order)
        for node in tree.iter_nodes()
    }
    resolved: dict[str, str] = {}
    root_signature = node_signature(tree.root)
    root_candidates = candidate_sets[root_signature]
    resolved[root_signature] = _best_supported_state(
        root_candidates,
        support_by_node[root_signature],
        dict.fromkeys(state_order, 1.0),
    )

    def visit(node, parent_state: str) -> None:
        for child in node.children:
            signature = node_signature(child)
            candidates = candidate_sets[signature]
            if parent_state in candidates:
                chosen = parent_state
            else:
                chosen = _best_supported_state(
                    candidates,
                    support_by_node[signature],
                    dict.fromkeys(state_order, 1.0),
                )
            resolved[signature] = chosen
            if not child.is_leaf():
                visit(child, chosen)

    visit(tree.root, resolved[root_signature])
    return resolved


def _transition_events(tree, resolved_states: dict[str, str]) -> list[TransitionEvent]:
    events: list[TransitionEvent] = []

    def visit(node) -> None:
        parent_signature = node_signature(node)
        parent_state = resolved_states[parent_signature]
        for child in node.children:
            child_signature = node_signature(child)
            child_state = resolved_states[child_signature]
            events.append(
                TransitionEvent(
                    parent_node=parent_signature,
                    child_node=child_signature,
                    source_state=parent_state,
                    target_state=child_state,
                    changed=parent_state != child_state,
                )
            )
            if not child.is_leaf():
                visit(child)

    visit(tree.root)
    return events


def _stationary_frequencies(
    states_by_taxon: dict[str, str], state_order: list[str]
) -> dict[str, float]:
    total = len(states_by_taxon)
    if total == 0:
        return dict.fromkeys(state_order, 0.0)
    return {
        state: float(
            format(
                sum(1 for value in states_by_taxon.values() if value == state) / total,
                ".15g",
            )
        )
        for state in state_order
    }


def _build_transition_count_matrix(
    state_order: list[str],
    er_events: list[TransitionEvent],
    *,
    model: str,
    state_ordering: str,
) -> dict[str, dict[str, float]]:
    counts = {
        source: {
            target: (
                1.0
                if target != source
                and _transition_allowed(
                    source, target, state_order, state_ordering=state_ordering
                )
                else 0.0
            )
            for target in state_order
            if target != source
        }
        for source in state_order
    }
    for event in er_events:
        if event.source_state == event.target_state:
            continue
        if not _transition_allowed(
            event.source_state,
            event.target_state,
            state_order,
            state_ordering=state_ordering,
        ):
            continue
        if model == "symmetric":
            counts[event.source_state][event.target_state] += 1.0
            if _transition_allowed(
                event.target_state,
                event.source_state,
                state_order,
                state_ordering=state_ordering,
            ):
                counts[event.target_state][event.source_state] += 1.0
        else:
            counts[event.source_state][event.target_state] += 1.0
    return counts


def _normal_interval(estimate: float, sample_size: float) -> tuple[float, float]:
    if sample_size <= 0.0:
        return estimate, estimate
    standard_error = math.sqrt(max(estimate * (1.0 - estimate), 0.0) / sample_size)
    lower = max(0.0, estimate - 1.96 * standard_error)
    upper = min(1.0, estimate + 1.96 * standard_error)
    return float(format(lower, ".15g")), float(format(upper, ".15g"))


def _fit_transition_matrix(
    model: str,
    state_order: list[str],
    stationary: dict[str, float],
    er_events: list[TransitionEvent],
    *,
    state_ordering: str,
) -> list[TransitionRateRow]:
    change_mass = 0.2
    stay_mass = 0.8
    if model == "equal-rates":
        return [
            TransitionRateRow(
                source_state=source,
                target_rates={
                    target: float(
                        format(
                            stay_mass
                            if source == target
                            else (
                                change_mass
                                / max(
                                    sum(
                                        1
                                        for candidate in state_order
                                        if candidate != source
                                        and _transition_allowed(
                                            source,
                                            candidate,
                                            state_order,
                                            state_ordering=state_ordering,
                                        )
                                    ),
                                    1,
                                )
                                if _transition_allowed(
                                    source,
                                    target,
                                    state_order,
                                    state_ordering=state_ordering,
                                )
                                else 0.0
                            ),
                            ".15g",
                        )
                    )
                    for target in state_order
                },
            )
            for source in state_order
        ]

    counts = _build_transition_count_matrix(
        state_order,
        er_events,
        model=model,
        state_ordering=state_ordering,
    )
    rows: list[TransitionRateRow] = []
    for source in state_order:
        off_total = sum(counts[source].values())
        rows.append(
            TransitionRateRow(
                source_state=source,
                target_rates={
                    target: (
                        float(format(stay_mass, ".15g"))
                        if source == target
                        else float(
                            format(
                                (
                                    change_mass * (counts[source][target] / off_total)
                                    if _transition_allowed(
                                        source,
                                        target,
                                        state_order,
                                        state_ordering=state_ordering,
                                    )
                                    and off_total > 0.0
                                    else 0.0
                                ),
                                ".15g",
                            )
                        )
                    )
                    for target in state_order
                },
            )
        )
    return rows


def _row_lookup(rows: list[TransitionRateRow]) -> dict[str, dict[str, float]]:
    return {row.source_state: row.target_rates for row in rows}


def _estimate_node_states(
    tree,
    candidate_sets: dict[str, list[str]],
    states_by_taxon: dict[str, str],
    state_order: list[str],
    matrix: list[TransitionRateRow],
    root_prior: dict[str, float],
    *,
    state_ordering: str,
) -> list[NodeStateEstimate]:
    transition_lookup = _row_lookup(matrix)
    support_by_node = {
        node_signature(node): _state_support(node, states_by_taxon, state_order)
        for node in tree.iter_nodes()
    }
    probabilities_by_node: dict[str, dict[str, float]] = {}
    resolved_states: dict[str, str] = {}

    root_signature = node_signature(tree.root)
    root_scores = {
        state: (support_by_node[root_signature].get(state, 0) + 1.0)
        * root_prior.get(state, 0.0)
        for state in (
            state_order
            if state_ordering == "ordered"
            else candidate_sets[root_signature]
        )
    }
    probabilities_by_node[root_signature] = _normalize_probabilities(root_scores)
    resolved_states[root_signature] = max(
        sorted(probabilities_by_node[root_signature]),
        key=lambda state: probabilities_by_node[root_signature][state],
    )

    def visit(node) -> None:
        parent_signature = node_signature(node)
        parent_probabilities = probabilities_by_node[parent_signature]
        for child in node.children:
            child_signature = node_signature(child)
            if child.is_leaf():
                leaf_state = states_by_taxon[child.name]
                probabilities_by_node[child_signature] = {leaf_state: 1.0}
                resolved_states[child_signature] = leaf_state
            else:
                scores: dict[str, float] = {}
                child_support = support_by_node[child_signature]
                for child_state in (
                    state_order
                    if state_ordering == "ordered"
                    else candidate_sets[child_signature]
                ):
                    transition_support = sum(
                        parent_probabilities[parent_state]
                        * transition_lookup[parent_state][child_state]
                        for parent_state in parent_probabilities
                    )
                    scores[child_state] = transition_support * (
                        child_support.get(child_state, 0) + 1.0
                    )
                probabilities_by_node[child_signature] = _normalize_probabilities(
                    scores
                )
                resolved_states[child_signature] = max(
                    sorted(probabilities_by_node[child_signature]),
                    key=lambda state: probabilities_by_node[child_signature][state],
                )
                visit(child)

    visit(tree.root)
    estimates: list[NodeStateEstimate] = []
    for node in tree.iter_nodes():
        signature = node_signature(node)
        state_probabilities = probabilities_by_node[signature]
        estimates.append(
            NodeStateEstimate(
                node=signature,
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=node_descendant_taxa(node),
                most_likely_state=resolved_states[signature],
                state_probabilities=state_probabilities,
                ambiguous=sum(
                    1
                    for probability in state_probabilities.values()
                    if probability > 0.0
                )
                > 1,
            )
        )
    return estimates


def _root_prior(
    model: str, stationary: dict[str, float], candidate_states: list[str]
) -> dict[str, float]:
    if model == "equal-rates":
        return dict.fromkeys(candidate_states, 1.0)
    return {state: stationary.get(state, 0.0) + 1e-9 for state in candidate_states}


def _pseudo_log_likelihood(
    estimates: list[NodeStateEstimate],
    events: list[TransitionEvent],
    model: TransitionModelReport,
) -> float:
    estimate_by_node = {estimate.node: estimate for estimate in estimates}
    transition_lookup = _row_lookup(model.transition_matrix)
    root_state = estimates[0].most_likely_state
    log_likelihood = math.log(
        max(model.root_state_probabilities.get(root_state, 1e-12), 1e-12)
    )
    for event in events:
        probability = transition_lookup[event.source_state][event.target_state]
        log_likelihood += math.log(max(probability, 1e-12))
        child_estimate = estimate_by_node[event.child_node]
        log_likelihood += math.log(
            max(
                child_estimate.state_probabilities.get(event.target_state, 1e-12), 1e-12
            )
        )
    return float(format(log_likelihood, ".15g"))


def _estimate_transition_rate_uncertainty(
    *,
    model: str,
    state_ordering: str,
    transition_matrix: list[TransitionRateRow],
    count_matrix: dict[str, dict[str, float]],
) -> TransitionRateUncertaintyReport:
    rows: list[TransitionRateUncertaintyRow] = []
    for row in transition_matrix:
        off_total = sum(count_matrix.get(row.source_state, {}).values())
        for target_state, estimate in row.target_rates.items():
            if row.source_state == target_state:
                lower, upper = estimate, estimate
                effective_count = off_total
            else:
                if off_total <= 0.0:
                    lower, upper = estimate, estimate
                else:
                    proportion = (
                        count_matrix[row.source_state][target_state] / off_total
                    )
                    lower_p, upper_p = _normal_interval(proportion, off_total)
                    lower = float(format(0.2 * lower_p, ".15g"))
                    upper = float(format(0.2 * upper_p, ".15g"))
                effective_count = count_matrix[row.source_state][target_state]
            rows.append(
                TransitionRateUncertaintyRow(
                    source_state=row.source_state,
                    target_state=target_state,
                    estimate=estimate,
                    lower_95_interval=lower,
                    upper_95_interval=upper,
                    effective_transition_count=effective_count,
                )
            )
    return TransitionRateUncertaintyReport(
        model=model,
        state_ordering=state_ordering,
        rows=rows,
    )


def _detect_sparse_state_instability(
    *,
    state_counts: dict[str, int],
    count_matrix: dict[str, dict[str, float]],
) -> SparseStateInstabilityReport:
    sparse_states = sorted(state for state, count in state_counts.items() if count < 2)
    zero_support_transitions = sorted(
        f"{source}->{target}"
        for source, targets in count_matrix.items()
        for target, value in targets.items()
        if value == 0.0
    )
    warning_count = int(bool(sparse_states)) + int(bool(zero_support_transitions))
    return SparseStateInstabilityReport(
        sparse_states=sparse_states,
        zero_support_transitions=zero_support_transitions,
        warning_count=warning_count,
        unstable=bool(sparse_states or zero_support_transitions),
    )


def _summarize_dominant_state_bias(
    state_counts: dict[str, int],
) -> DominantStateBiasReport:
    if not state_counts:
        return DominantStateBiasReport(
            dominant_states=[],
            dominant_fraction=0.0,
            biased=False,
            message=None,
        )
    total = sum(state_counts.values())
    max_count = max(state_counts.values())
    dominant_states = sorted(
        state for state, count in state_counts.items() if count == max_count
    )
    dominant_fraction = float(format(max_count / max(total, 1), ".15g"))
    biased = dominant_fraction >= 0.8
    return DominantStateBiasReport(
        dominant_states=dominant_states,
        dominant_fraction=dominant_fraction,
        biased=biased,
        message=(
            "one state dominates most taxa and may compress minority-state transition evidence"
            if biased
            else None
        ),
    )


def assess_geographic_state_analysis_readiness(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> GeographicAnalysisReadinessReport:
    """Decide whether one geographic discrete-state analysis is credible enough to run."""
    tree_validation = validate_tree_path(tree_path)
    coding = validate_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    imbalance = detect_state_imbalance_problems(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    dominant_state_bias = _summarize_dominant_state_bias(imbalance.state_counts)
    blockers: list[str] = []
    warnings: list[str] = []

    if not tree_validation.syntax_valid or not tree_validation.biologically_safe:
        blockers.append(
            "tree validation failed and the geographic analysis is not safe to interpret"
        )
    if not tree_validation.rooted:
        blockers.append("geographic ancestral-state analysis requires a rooted tree")
    if not coding.valid:
        blockers.append(
            "discrete geographic states contain unsupported labels or coding patterns"
        )
    if any(warning.code == "single-state-dataset" for warning in imbalance.warnings):
        blockers.append(
            "geographic analysis requires at least two observed states after matching taxa to the tree"
        )
    rare_state_count = sum(1 for count in imbalance.state_counts.values() if count < 2)
    if imbalance.state_counts and rare_state_count == len(imbalance.state_counts):
        blockers.append(
            "one or more geographic states are too sparse to estimate transitions credibly"
        )
    if dominant_state_bias.biased:
        blockers.append(
            "observed geographic states are dominated by one state and the sampling is too biased for credible transition inference"
        )

    warnings.extend(tree_validation.warnings)
    warnings.extend(warning.message for warning in imbalance.warnings)
    if (
        dominant_state_bias.message is not None
        and dominant_state_bias.message not in warnings
    ):
        warnings.append(dominant_state_bias.message)

    return GeographicAnalysisReadinessReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=coding.taxon_column,
        trait=trait,
        valid=not blockers,
        blockers=blockers,
        warnings=warnings,
        state_ordering=state_ordering,
        ordered_states=list(ordered_states or []),
        coding_validation=coding,
        imbalance=imbalance,
        dominant_state_bias=dominant_state_bias,
        tree_validation_decision=tree_validation.validity_decision,
    )


def _estimate_transition_support_rows(
    *,
    estimates: list[NodeStateEstimate],
    events: list[TransitionEvent],
    transition_matrix: list[TransitionRateRow],
) -> list[TransitionSupportRow]:
    estimate_by_node = {estimate.node: estimate for estimate in estimates}
    transition_lookup = _row_lookup(transition_matrix)
    rows: list[TransitionSupportRow] = []
    for event in events:
        parent_probabilities = estimate_by_node[event.parent_node].state_probabilities
        child_probabilities = estimate_by_node[event.child_node].state_probabilities
        scores = {
            (parent_state, child_state): (
                parent_probability
                * transition_lookup[parent_state][child_state]
                * child_probabilities.get(child_state, 0.0)
            )
            for parent_state, parent_probability in parent_probabilities.items()
            for child_state in child_probabilities
        }
        total_score = sum(scores.values())
        inferred_key = (event.source_state, event.target_state)
        support = (
            0.0 if total_score <= 0.0 else scores.get(inferred_key, 0.0) / total_score
        )
        rows.append(
            TransitionSupportRow(
                parent_node=event.parent_node,
                child_node=event.child_node,
                inferred_transition=f"{event.source_state}->{event.target_state}",
                support=float(format(support, ".15g")),
                strongly_supported=event.changed and support >= 0.8,
            )
        )
    return rows


def _quantile(sorted_values: list[float], fraction: float) -> float:
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return float(format(sorted_values[0], ".15g"))
    index = max(
        0, min(len(sorted_values) - 1, int(round(fraction * (len(sorted_values) - 1))))
    )
    return float(format(sorted_values[index], ".15g"))


def _build_discrete_evolution_narrative(
    report: DiscreteStateEvolutionReport,
    *,
    comparison: DiscreteModelComparisonReport | None = None,
) -> DiscreteEvolutionNarrative:
    root_state = next(
        estimate.most_likely_state
        for estimate in report.estimates
        if estimate.node == node_signature(load_tree(report.tree_path).root)
    )
    summary = (
        f"reconstructed {report.trait} across {report.taxon_count} taxa under the {report.model} model"
        f" with root state '{root_state}' and {report.transition_summary.transition_count} inferred branch transitions"
    )
    strong = report.transition_summary.strongly_supported_transition_count
    transition_summary = (
        f"{strong} of {report.transition_summary.transition_count} changed branches remain strongly supported after"
        " weighting by parent-child state probabilities"
    )
    caveats = list(report.warnings)
    if comparison is not None and comparison.sensitive_region_count > 0:
        caveats.append(
            f"{comparison.sensitive_region_count} internal regions change reconstructed state between {comparison.left_model}"
            f" and {comparison.right_model}"
        )
    caveats.extend(
        [
            "deterministic node probabilities are an approximation and do not replace full Bayesian or likelihood-marginal ancestral mapping",
            "transition uncertainty intervals summarize effective event counts in the fitted deterministic model and are not posterior credible intervals",
        ]
    )
    return DiscreteEvolutionNarrative(
        summary=summary,
        transition_summary=transition_summary,
        interpretation_boundary=(
            "treat these outputs as computational evidence about state histories, not as direct proof of dispersal timing, mechanism, or causal biogeography"
        ),
        caveats=caveats,
    )


def _build_model_sensitive_regions(
    differences: list[NodeStateDifference],
) -> list[ModelSensitiveRegionRow]:
    rows: list[ModelSensitiveRegionRow] = []
    for difference in differences:
        if not difference.differs:
            continue
        left_probability = difference.left_probabilities.get(difference.left_state, 0.0)
        right_probability = difference.right_probabilities.get(
            difference.right_state, 0.0
        )
        rows.append(
            ModelSensitiveRegionRow(
                node=difference.node,
                descendant_taxa=difference.descendant_taxa,
                left_state=difference.left_state,
                right_state=difference.right_state,
                sensitivity_score=float(
                    format(abs(left_probability - right_probability), ".15g")
                ),
            )
        )
    return sorted(rows, key=lambda row: (-row.sensitivity_score, row.node))


def _sample_state(probabilities: dict[str, float], rng: random.Random) -> str:
    threshold = rng.random()
    cumulative = 0.0
    ordered_items = sorted(probabilities.items())
    for state, probability in ordered_items:
        cumulative += probability
        if threshold <= cumulative:
            return state
    return ordered_items[-1][0]


def _sample_index(weights: numpy.ndarray, rng: random.Random) -> int:
    total = float(weights.sum())
    if total <= 0.0:
        raise AncestralReconstructionError(
            "cannot sample a discrete state from zero-probability weights"
        )
    threshold = rng.random() * total
    cumulative = 0.0
    for index, weight in enumerate(weights):
        cumulative += float(weight)
        if threshold <= cumulative:
            return index
    return len(weights) - 1


def _sample_transition_target(
    source_state: str,
    transition_lookup: dict[str, dict[str, float]],
    rng: random.Random,
) -> str:
    off_diagonal = {
        target_state: probability
        for target_state, probability in transition_lookup[source_state].items()
        if target_state != source_state and probability > 0.0
    }
    if not off_diagonal:
        return source_state
    total = sum(off_diagonal.values())
    normalized = {
        state: float(format(probability / total, ".15g"))
        for state, probability in off_diagonal.items()
    }
    return _sample_state(normalized, rng)


def _postorder_discrete_partials(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
) -> dict[str, numpy.ndarray]:
    state_index = {state: index for index, state in enumerate(state_order)}
    transition_cache: dict[float, numpy.ndarray] = {}
    partials: dict[str, numpy.ndarray] = {}

    def transition(branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(branch_length)
        if cached is None:
            cached = _transition_probability_matrix(rate_matrix, branch_length)
            transition_cache[branch_length] = cached
        return cached

    def visit(node) -> numpy.ndarray:
        signature = node_signature(node)
        if node.is_leaf():
            partial = numpy.zeros(len(state_order), dtype=float)
            partial[state_index[states_by_taxon[node.name]]] = 1.0
            partials[signature] = partial
            return partial
        partial = numpy.ones(len(state_order), dtype=float)
        for child in node.children:
            partial *= transition(_branch_length(child)) @ visit(child)
        partials[signature] = partial
        return partial

    visit(tree.root)
    return partials


def _sample_conditional_node_states(
    tree,
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
    partials: dict[str, numpy.ndarray],
    rng: random.Random,
) -> dict[str, str]:
    state_index = {state: index for index, state in enumerate(state_order)}
    transition_cache: dict[float, numpy.ndarray] = {}
    node_states: dict[str, str] = {}

    def transition(branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(branch_length)
        if cached is None:
            cached = _transition_probability_matrix(rate_matrix, branch_length)
            transition_cache[branch_length] = cached
        return cached

    root_signature = node_signature(tree.root)
    root_weights = root_prior * partials[root_signature]
    node_states[root_signature] = state_order[_sample_index(root_weights, rng)]

    def visit(node) -> None:
        parent_signature = node_signature(node)
        parent_index = state_index[node_states[parent_signature]]
        for child in node.children:
            child_signature = node_signature(child)
            child_weights = (
                partials[child_signature]
                if child.is_leaf()
                else transition(_branch_length(child))[parent_index, :]
                * partials[child_signature]
            )
            node_states[child_signature] = state_order[_sample_index(child_weights, rng)]
            if not child.is_leaf():
                visit(child)

    visit(tree.root)
    return node_states


def _simulate_ctmc_branch_path(
    *,
    branch_index: int,
    parent_node: str,
    child_node: str,
    branch_length: float,
    start_state: str,
    end_state: str,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    rng: random.Random,
    max_attempts: int = 2000,
) -> tuple[list[StochasticMapTransitionEvent], list[StochasticMapStateSegment], int]:
    if branch_length <= 0.0:
        if start_state != end_state:
            raise AncestralReconstructionError(
                "zero-length branch cannot support different start and end states"
            )
        return (
            [],
            [
                StochasticMapStateSegment(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    state=start_state,
                    start_time_fraction=0.0,
                    end_time_fraction=1.0,
                    duration=0.0,
                )
            ],
            0,
        )

    state_index = {state: index for index, state in enumerate(state_order)}

    def choose_target(current_index: int) -> int:
        row = rate_matrix[current_index, :]
        weights = numpy.array(
            [
                0.0 if index == current_index else float(rate)
                for index, rate in enumerate(row)
            ],
            dtype=float,
        )
        return _sample_index(weights, rng)

    for attempt_index in range(1, max_attempts + 1):
        current_state = start_state
        current_index = state_index[current_state]
        elapsed = 0.0
        events: list[StochasticMapTransitionEvent] = []
        while elapsed < branch_length:
            exit_rate = float(-rate_matrix[current_index, current_index])
            if exit_rate <= 0.0:
                break
            elapsed += rng.expovariate(exit_rate)  # nosec B311
            if elapsed >= branch_length:
                break
            next_index = choose_target(current_index)
            next_state = state_order[next_index]
            if next_state == current_state:
                continue
            events.append(
                StochasticMapTransitionEvent(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    source_state=current_state,
                    target_state=next_state,
                    event_time_fraction=float(format(elapsed / branch_length, ".15g")),
                )
            )
            current_state = next_state
            current_index = next_index
        if current_state != end_state:
            continue
        segments: list[StochasticMapStateSegment] = []
        segment_state = start_state
        segment_start = 0.0
        for event in events:
            segment_end = min(max(event.event_time_fraction, segment_start), 1.0)
            segments.append(
                StochasticMapStateSegment(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    state=segment_state,
                    start_time_fraction=float(format(segment_start, ".15g")),
                    end_time_fraction=float(format(segment_end, ".15g")),
                    duration=float(
                        format(branch_length * (segment_end - segment_start), ".15g")
                    ),
                )
            )
            segment_state = event.target_state
            segment_start = segment_end
        segments.append(
            StochasticMapStateSegment(
                branch_index=branch_index,
                parent_node=parent_node,
                child_node=child_node,
                state=segment_state,
                start_time_fraction=float(format(segment_start, ".15g")),
                end_time_fraction=1.0,
                duration=float(
                    format(branch_length * max(1.0 - segment_start, 0.0), ".15g")
                ),
            )
        )
        return events, segments, attempt_index
    raise AncestralReconstructionError(
        "failed to sample a branch history consistent with the conditioned endpoint states"
    )


def _summarize_stochastic_map_replicates(
    replicates: list[StochasticMapReplicate],
    *,
    simulation_failure_count: int,
    expected_transitions: list[str] | None = None,
) -> StochasticMapSummaryReport:
    total_counts = sorted(
        float(replicate.total_transition_count) for replicate in replicates
    )
    transition_names = sorted(
        set(expected_transitions or [])
        | {
            transition
            for replicate in replicates
            for transition in replicate.transition_counts
        }
    )
    rows: list[StochasticMapSummaryRow] = []
    for transition in transition_names:
        values = [
            replicate.transition_counts.get(transition, 0) for replicate in replicates
        ]
        sorted_values = sorted(float(value) for value in values)
        rows.append(
            StochasticMapSummaryRow(
                transition=transition,
                mean_count=float(format(sum(values) / max(len(values), 1), ".15g")),
                lower_95_interval=_quantile(sorted_values, 0.025),
                upper_95_interval=_quantile(sorted_values, 0.975),
                minimum_count=min(values, default=0),
                maximum_count=max(values, default=0),
                presence_fraction=float(
                    format(
                        sum(1 for value in values if value > 0) / max(len(values), 1),
                        ".15g",
                    )
                ),
            )
        )
    state_names = sorted(
        {
            state
            for replicate in replicates
            for state in replicate.state_time_totals
        }
    )
    state_time_rows: list[StochasticMapStateTimeRow] = []
    for state in state_names:
        values = [
            float(replicate.state_time_totals.get(state, 0.0))
            for replicate in replicates
        ]
        sorted_values = sorted(values)
        state_time_rows.append(
            StochasticMapStateTimeRow(
                state=state,
                mean_time=float(format(sum(values) / max(len(values), 1), ".15g")),
                lower_95_interval=_quantile(sorted_values, 0.025),
                upper_95_interval=_quantile(sorted_values, 0.975),
                minimum_time=min(values, default=0.0),
                maximum_time=max(values, default=0.0),
            )
        )
    branch_lookup: dict[tuple[int, str, str, float], list[dict[str, float]]] = {}
    for replicate in replicates:
        for history in replicate.branch_histories:
            key = (
                history.branch_index,
                history.parent_node,
                history.child_node,
                float(history.branch_length),
            )
            branch_lookup.setdefault(key, [])
            state_times = {
                state: 0.0
                for state in state_names
            }
            for segment in history.segments:
                state_times[segment.state] = state_times.get(segment.state, 0.0) + float(
                    segment.duration
                )
            branch_lookup[key].append(state_times)
    branch_occupancy_rows: list[StochasticMapBranchOccupancyRow] = []
    for (
        branch_index,
        parent_node,
        child_node,
        branch_length,
    ), state_times in sorted(
        branch_lookup.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
            item[0][2],
        ),
    ):
        for state in state_names:
            values = [
                replicate_state_times.get(state, 0.0)
                for replicate_state_times in state_times
            ]
            sorted_values = sorted(float(value) for value in values)
            mean_time = float(format(sum(values) / max(len(values), 1), ".15g"))
            mean_fraction = 0.0
            if branch_length > 0.0:
                mean_fraction = float(format(mean_time / branch_length, ".15g"))
            branch_occupancy_rows.append(
                StochasticMapBranchOccupancyRow(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    state=state,
                    branch_length=branch_length,
                    mean_time=mean_time,
                    lower_95_interval=_quantile(sorted_values, 0.025),
                    upper_95_interval=_quantile(sorted_values, 0.975),
                    minimum_time=min(values, default=0.0),
                    maximum_time=max(values, default=0.0),
                    mean_fraction=mean_fraction,
                    presence_fraction=float(
                        format(
                            sum(1 for value in values if value > 0.0)
                            / max(len(values), 1),
                            ".15g",
                        )
                    ),
                )
            )
    warnings: list[str] = []
    if simulation_failure_count > 0:
        warnings.append(
            "one or more stochastic-map replicates failed to sample a branch history consistent with the conditioned endpoint states"
        )
    return StochasticMapSummaryReport(
        replicate_count=len(replicates),
        mean_total_transition_count=float(
            format(sum(total_counts) / max(len(total_counts), 1), ".15g")
        ),
        lower_95_total_transition_count=_quantile(total_counts, 0.025),
        upper_95_total_transition_count=_quantile(total_counts, 0.975),
        rows=rows,
        state_time_rows=state_time_rows,
        branch_occupancy_rows=branch_occupancy_rows,
        simulation_failure_count=simulation_failure_count,
        warnings=warnings,
    )


def _stochastic_map_warning_union(*warning_lists: list[str]) -> list[str]:
    merged: list[str] = []
    for warnings in warning_lists:
        for warning in warnings:
            if warning not in merged:
                merged.append(warning)
    return merged


def _stochastic_map_fit_audit(
    fit_report: DiscreteMkFitReport,
) -> StochasticMapModelFitAudit:
    baseline = fit_report.baseline_comparison
    diagnostics = fit_report.optimizer_diagnostics
    return StochasticMapModelFitAudit(
        state_order=list(fit_report.state_order),
        allowed_transitions=[
            f"{row.source_state}->{row.target_state}"
            for row in fit_report.transition_rate_rows
            if row.transition_allowed and row.source_state != row.target_state
        ],
        parameter_count=fit_report.parameter_count,
        log_likelihood=fit_report.log_likelihood,
        aic=fit_report.aic,
        aicc=fit_report.aicc,
        overparameterized=fit_report.overparameterized,
        optimizer_converged=diagnostics.converged,
        optimizer_iteration_count=diagnostics.iteration_count,
        optimizer_function_evaluation_count=diagnostics.function_evaluation_count,
        optimizer_hit_lower_parameter_bound=diagnostics.hit_lower_parameter_bound,
        optimizer_hit_upper_parameter_bound=diagnostics.hit_upper_parameter_bound,
        baseline_model=(None if baseline is None else baseline.baseline_model),
        baseline_aic=(None if baseline is None else baseline.baseline_aic),
        baseline_delta_aic=(None if baseline is None else baseline.delta_aic),
        preferred_model_by_aic=(
            None if baseline is None else baseline.preferred_model_by_aic
        ),
        warnings=list(fit_report.input_audit.warnings),
    )


def _rate_matrix_from_transition_rate_rows(
    *,
    state_order: list[str],
    transition_rate_rows,
) -> numpy.ndarray:
    state_index = {state: index for index, state in enumerate(state_order)}
    rate_matrix = numpy.zeros((len(state_order), len(state_order)), dtype=float)
    for row in transition_rate_rows:
        if not row.transition_allowed:
            continue
        if row.source_state == row.target_state:
            continue
        rate_matrix[
            state_index[row.source_state],
            state_index[row.target_state],
        ] = float(row.rate)
    numpy.fill_diagonal(rate_matrix, -numpy.sum(rate_matrix, axis=1))
    return rate_matrix


def _simulate_stochastic_maps_from_components(
    *,
    dataset,
    model: str,
    state_order: list[str],
    state_ordering: str,
    ordered_states: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
    replicates: int,
    seed: int,
    fit_audit: StochasticMapModelFitAudit,
) -> StochasticMapCollectionReport:
    partials = _postorder_discrete_partials(
        dataset.tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
    )
    branch_rows = [
        (
            index,
            node_signature(node.parent),
            node_signature(node),
            node,
        )
        for index, node in enumerate(dataset.tree.iter_nodes())
        if node.parent is not None
    ]
    randomizer = random.Random(seed)  # nosec B311
    maps: list[StochasticMapReplicate] = []
    failures: list[StochasticMapSimulationFailure] = []
    for replicate_index in range(replicates):
        node_states = _sample_conditional_node_states(
            dataset.tree,
            state_order=state_order,
            rate_matrix=rate_matrix,
            root_prior=root_prior,
            partials=partials,
            rng=randomizer,
        )
        root_state = node_states[node_signature(dataset.tree.root)]
        branch_histories: list[StochasticMapBranchHistory] = []
        transition_counts: dict[str, int] = {}
        state_time_totals = {state: 0.0 for state in state_order}
        total_transition_count = 0
        for branch_index, parent_node, child_node, child in branch_rows:
            parent_state = node_states[parent_node]
            child_state = node_states[child_node]
            branch_length = _branch_length(child)
            try:
                events, segments, attempt_count = _simulate_ctmc_branch_path(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    branch_length=branch_length,
                    start_state=parent_state,
                    end_state=child_state,
                    state_order=state_order,
                    rate_matrix=rate_matrix,
                    rng=randomizer,
                )
            except AncestralReconstructionError:
                failures.append(
                    StochasticMapSimulationFailure(
                        replicate_index=replicate_index,
                        branch_index=branch_index,
                        parent_node=parent_node,
                        child_node=child_node,
                        source_state=parent_state,
                        target_state=child_state,
                        branch_length=branch_length,
                        attempt_count=2000,
                        reason=(
                            "failed to sample a branch history consistent with the conditioned endpoint states"
                        ),
                    )
                )
                branch_histories = []
                transition_counts = {}
                state_time_totals = {state: 0.0 for state in state_order}
                total_transition_count = 0
                break
            for event in events:
                transition = f"{event.source_state}->{event.target_state}"
                transition_counts[transition] = transition_counts.get(transition, 0) + 1
                total_transition_count += 1
            for segment in segments:
                state_time_totals[segment.state] = float(
                    format(
                        state_time_totals.get(segment.state, 0.0) + segment.duration,
                        ".15g",
                    )
                )
            branch_histories.append(
                StochasticMapBranchHistory(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    branch_length=branch_length,
                    start_state=parent_state,
                    end_state=child_state,
                    event_count=len(events),
                    events=events,
                    segments=segments,
                )
            )
        if not branch_histories:
            continue
        maps.append(
            StochasticMapReplicate(
                replicate_index=replicate_index,
                root_state=root_state,
                total_transition_count=total_transition_count,
                transition_counts=dict(sorted(transition_counts.items())),
                state_time_totals=dict(sorted(state_time_totals.items())),
                branch_histories=branch_histories,
            )
        )
    if not maps:
        raise AncestralReconstructionError(
            "stochastic character mapping failed for every requested replicate"
        )
    summary = _summarize_stochastic_map_replicates(
        maps,
        simulation_failure_count=len(failures),
        expected_transitions=fit_audit.allowed_transitions,
    )
    return StochasticMapCollectionReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        model=model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        replicates=replicates,
        seed=seed,
        conditioned_on_node_estimates=False,
        fit_audit=fit_audit,
        warnings=_stochastic_map_warning_union(fit_audit.warnings, summary.warnings),
        maps=maps,
        failures=failures,
        summary=summary,
    )


def validate_discrete_state_coding(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> StateCodingValidationReport:
    """Detect impossible or unsupported discrete-state labels."""
    tree = load_tree(tree_path)
    table = (
        load_tsv_summary(traits_path)
        if taxon_column is None
        else load_taxon_table(traits_path, taxon_column=taxon_column)
    )
    if trait not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{trait}'"
        )
    tree_taxa = set(tree.tip_names)
    if state_ordering not in {"unordered", "ordered"}:
        raise ValueError(f"unsupported state ordering: {state_ordering}")
    ordered = list(ordered_states or [])
    if ordered and len(set(ordered)) != len(ordered):
        raise AncestralReconstructionError(
            "ordered state vocabulary contains duplicate labels"
        )
    allowed = list(allowed_states or ordered or [])
    allowed_set = set(allowed)
    issues: list[StateCodingIssue] = []
    usable_taxa: list[str] = []
    observed_states: set[str] = set()
    for row in table.rows:
        taxon = row[table.taxon_column]
        if taxon not in tree_taxa:
            continue
        raw_state = row[trait].strip()
        if not raw_state:
            continue
        invalid_delimiter = next(
            (
                token
                for token in _DEFAULT_ALLOWED_STATE_PATTERN_BLOCKLIST
                if token in raw_state
            ),
            None,
        )
        if invalid_delimiter is not None:
            issues.append(
                StateCodingIssue(
                    taxon=taxon,
                    raw_state=raw_state,
                    code="unsupported-state-delimiter",
                    message=f"state label contains reserved delimiter '{invalid_delimiter}'",
                )
            )
            continue
        if allowed and raw_state not in allowed_set:
            issues.append(
                StateCodingIssue(
                    taxon=taxon,
                    raw_state=raw_state,
                    code="unordered-state-vocabulary"
                    if state_ordering == "ordered" and ordered
                    else "unsupported-state-label",
                    message=(
                        "state label is not present in the declared ordered vocabulary"
                        if state_ordering == "ordered" and ordered
                        else "state label is not present in the allowed state vocabulary"
                    ),
                )
            )
            continue
        observed_states.add(raw_state)
        usable_taxa.append(taxon)
    if state_ordering == "ordered" and ordered:
        missing_from_order = sorted(observed_states - set(ordered))
        issues.extend(
            StateCodingIssue(
                taxon="",
                raw_state=state,
                code="unordered-state-vocabulary",
                message="observed state is missing from the declared ordered vocabulary",
            )
            for state in missing_from_order
        )
    return StateCodingValidationReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        allowed_states=allowed,
        state_ordering=state_ordering,
        ordered_states=ordered,
        valid=not issues,
        issues=issues,
        observed_states=sorted(observed_states),
        usable_taxa=sorted(usable_taxa),
    )


def audit_discrete_state_coding(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    coding_map: dict[str, str] | None = None,
) -> StateCodingAuditReport:
    """Show how raw metadata states become model states or get excluded."""
    tree = load_tree(tree_path)
    table = (
        load_tsv_summary(traits_path)
        if taxon_column is None
        else load_taxon_table(traits_path, taxon_column=taxon_column)
    )
    if trait not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{trait}'"
        )
    tree_taxa = set(tree.tip_names)
    mapping = dict(coding_map or {})
    ordered = list(ordered_states or [])
    allowed = list(allowed_states or ordered or [])
    allowed_set = set(allowed)
    rows: list[StateCodingAuditRow] = []
    for row in table.rows:
        taxon = row[table.taxon_column]
        raw_state = row[trait].strip()
        in_tree = taxon in tree_taxa
        if not raw_state:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=None,
                    in_tree=in_tree,
                    included=False,
                    issue_code="missing-state",
                    note="state is blank and cannot be used",
                )
            )
            continue
        normalized_state = mapping.get(raw_state, raw_state)
        invalid_delimiter = next(
            (
                token
                for token in _DEFAULT_ALLOWED_STATE_PATTERN_BLOCKLIST
                if token in normalized_state
            ),
            None,
        )
        if not in_tree:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=normalized_state,
                    in_tree=False,
                    included=False,
                    issue_code="taxon-not-in-tree",
                    note="taxon does not overlap the tree and is excluded before state modeling",
                )
            )
            continue
        if invalid_delimiter is not None:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=normalized_state,
                    in_tree=True,
                    included=False,
                    issue_code="unsupported-state-delimiter",
                    note=f"normalized state contains reserved delimiter '{invalid_delimiter}'",
                )
            )
            continue
        if allowed and normalized_state not in allowed_set:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=normalized_state,
                    in_tree=True,
                    included=False,
                    issue_code="unsupported-state-label"
                    if state_ordering != "ordered"
                    else "unordered-state-vocabulary",
                    note=(
                        "normalized state is absent from the declared ordered vocabulary"
                        if state_ordering == "ordered" and ordered
                        else "normalized state is absent from the allowed state vocabulary"
                    ),
                )
            )
            continue
        rows.append(
            StateCodingAuditRow(
                taxon=taxon,
                raw_state=raw_state,
                normalized_state=normalized_state,
                in_tree=True,
                included=True,
                issue_code=None,
                note="state is retained for discrete-state modeling",
            )
        )
    included_row_count = sum(1 for row in rows if row.included)
    return StateCodingAuditReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        state_ordering=state_ordering,
        ordered_states=ordered,
        coding_map=mapping,
        row_count=len(rows),
        included_row_count=included_row_count,
        excluded_row_count=len(rows) - included_row_count,
        rows=rows,
    )


def detect_state_imbalance_problems(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> StateImbalanceReport:
    """Flag rare or degenerate discrete-state inputs over the tree overlap."""
    tree = load_tree(tree_path)
    table = (
        load_tsv_summary(traits_path)
        if taxon_column is None
        else load_taxon_table(traits_path, taxon_column=taxon_column)
    )
    if trait not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{trait}'"
        )
    tree_taxa = set(tree.tip_names)
    state_counts: dict[str, int] = {}
    for row in table.rows:
        taxon = row[table.taxon_column]
        state = row[trait].strip()
        if taxon in tree_taxa and state:
            state_counts[state] = state_counts.get(state, 0) + 1
    observed_states = sorted(state_counts)
    warnings: list[StateImbalanceWarning] = []
    rare_states = sorted(state for state, count in state_counts.items() if count < 2)
    if len(observed_states) < 2:
        warnings.append(
            StateImbalanceWarning(
                code="single-state-dataset",
                message="only one observed state remains after pruning to usable tree taxa",
                affected_states=observed_states,
            )
        )
    if rare_states:
        warnings.append(
            StateImbalanceWarning(
                code="rare-states",
                message="one or more states are represented by fewer than two taxa",
                affected_states=rare_states,
            )
        )
    dominant_fraction = (
        max(state_counts.values()) / max(sum(state_counts.values()), 1)
        if state_counts
        else 0.0
    )
    if dominant_fraction >= 0.8 and observed_states:
        dominant_states = [
            state
            for state, count in state_counts.items()
            if count == max(state_counts.values())
        ]
        warnings.append(
            StateImbalanceWarning(
                code="dominant-state-skew",
                message="one state dominates most observed taxa and may overwhelm transition inference",
                affected_states=sorted(dominant_states),
            )
        )
    return StateImbalanceReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        taxon_count=sum(state_counts.values()),
        observed_states=observed_states,
        state_counts=state_counts,
        warnings=warnings,
    )


def run_discrete_state_transition_model(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteStateEvolutionReport:
    """Run a deterministic discrete-state evolution workflow on one tree and trait."""
    if model not in {"equal-rates", "symmetric", "all-rates-different"}:
        raise ValueError(f"unsupported discrete-state model: {model}")
    coding = validate_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    if not coding.valid:
        raise AncestralReconstructionError(
            "discrete-state evolution input contains unsupported state labels"
        )
    imbalance = detect_state_imbalance_problems(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    if any(warning.code == "single-state-dataset" for warning in imbalance.warnings):
        raise AncestralReconstructionError(
            "discrete-state evolution requires at least two observed states"
        )
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    state_order = _resolve_state_order(
        dataset.observed_states,
        allowed_states=allowed_states,
        ordered_states=ordered_states,
        state_ordering=state_ordering,
    )
    candidate_sets = _fitch_candidate_sets(dataset.tree, dataset.states_by_taxon)
    stationary = _stationary_frequencies(dataset.states_by_taxon, state_order)
    er_resolved = _resolve_er_states(
        dataset.tree, candidate_sets, dataset.states_by_taxon, state_order
    )
    er_events = _transition_events(dataset.tree, er_resolved)
    count_matrix = _build_transition_count_matrix(
        state_order,
        er_events,
        model=model,
        state_ordering=state_ordering,
    )
    matrix = _fit_transition_matrix(
        model, state_order, stationary, er_events, state_ordering=state_ordering
    )
    root_prior = _root_prior(
        model, stationary, candidate_sets[node_signature(dataset.tree.root)]
    )
    estimates = _estimate_node_states(
        dataset.tree,
        candidate_sets,
        dataset.states_by_taxon,
        state_order,
        matrix,
        root_prior,
        state_ordering=state_ordering,
    )
    resolved_states = {
        estimate.node: estimate.most_likely_state for estimate in estimates
    }
    events = _transition_events(dataset.tree, resolved_states)
    transition_counts: dict[str, int] = {}
    for event in events:
        key = f"{event.source_state}->{event.target_state}"
        transition_counts[key] = transition_counts.get(key, 0) + 1
    support_rows = _estimate_transition_support_rows(
        estimates=estimates,
        events=events,
        transition_matrix=matrix,
    )
    branch_count = len(events)
    transition_count = sum(1 for event in events if event.changed)
    strongly_supported_transition_count = sum(
        1 for row in support_rows if row.strongly_supported
    )
    strongly_supported_transition_counts: dict[str, int] = {}
    for row in support_rows:
        if row.strongly_supported:
            strongly_supported_transition_counts[row.inferred_transition] = (
                strongly_supported_transition_counts.get(row.inferred_transition, 0) + 1
            )
    uncertainty = _estimate_transition_rate_uncertainty(
        model=model,
        state_ordering=state_ordering,
        transition_matrix=matrix,
        count_matrix=count_matrix,
    )
    instability = _detect_sparse_state_instability(
        state_counts=dataset.state_counts,
        count_matrix=count_matrix,
    )
    dominant_state_bias = _summarize_dominant_state_bias(dataset.state_counts)
    transition_model = TransitionModelReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=model,
        likelihood_method="deterministic-node-probability",
        state_ordering=state_ordering,
        ordered_states=state_order if state_ordering == "ordered" else [],
        state_order=state_order,
        parameter_count=(
            1
            if model == "equal-rates"
            else (
                len(state_order) * max(len(state_order) - 1, 0) // 2
                if model == "symmetric"
                else len(state_order) * max(len(state_order) - 1, 0)
            )
        ),
        pseudo_log_likelihood=0.0,
        aic=0.0,
        stationary_frequencies=stationary,
        transition_matrix=matrix,
        uncertainty=uncertainty,
        root_state_probabilities=_normalize_probabilities(root_prior),
    )
    transition_model.pseudo_log_likelihood = _pseudo_log_likelihood(
        estimates, events, transition_model
    )
    transition_model.aic = float(
        format(
            2.0 * transition_model.parameter_count
            - 2.0 * transition_model.pseudo_log_likelihood,
            ".15g",
        )
    )
    summary = TransitionSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=model,
        branch_count=branch_count,
        transition_count=transition_count,
        strongly_supported_transition_count=strongly_supported_transition_count,
        transition_counts=dict(sorted(transition_counts.items())),
        strongly_supported_transition_counts=dict(
            sorted(strongly_supported_transition_counts.items())
        ),
        support_rows=support_rows,
        events=events,
    )
    warnings = list(dataset.warnings)
    warnings.extend(warning.message for warning in imbalance.warnings)
    if instability.unstable:
        warnings.append(
            "sparse-state transition estimates may be unstable for one or more source-target paths"
        )
    if dominant_state_bias.biased and dominant_state_bias.message is not None:
        warnings.append(dominant_state_bias.message)
    return DiscreteStateEvolutionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=model,
        likelihood_method="deterministic-node-probability",
        state_ordering=state_ordering,
        ordered_states=state_order if state_ordering == "ordered" else [],
        analysis_tree_newick=dumps_newick(dataset.tree),
        taxon_count=len(dataset.taxa),
        observed_states=state_order,
        state_counts=dataset.state_counts,
        coding_validation=coding,
        imbalance=imbalance,
        instability=instability,
        dominant_state_bias=dominant_state_bias,
        transition_model=transition_model,
        estimates=estimates,
        transition_summary=summary,
        warnings=warnings,
    )


def estimate_ancestral_geographic_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteStateEvolutionReport:
    """Estimate ancestral geographic states over a rooted tree."""
    readiness = assess_geographic_state_analysis_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    if not readiness.valid:
        raise AncestralReconstructionError(
            "geographic state analysis is inappropriate: "
            + "; ".join(readiness.blockers)
        )
    return run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )


def compare_discrete_state_models(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    left_model: str = "equal-rates",
    right_model: str = "all-rates-different",
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteModelComparisonReport:
    """Compare discrete-state reconstructions across two supported models."""
    left = run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=left_model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    right = run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=right_model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    right_by_node = {estimate.node: estimate for estimate in right.estimates}
    differences: list[NodeStateDifference] = []
    for left_estimate in left.estimates:
        right_estimate = right_by_node[left_estimate.node]
        differences.append(
            NodeStateDifference(
                node=left_estimate.node,
                descendant_taxa=left_estimate.descendant_taxa,
                left_state=left_estimate.most_likely_state,
                right_state=right_estimate.most_likely_state,
                differs=left_estimate.most_likely_state
                != right_estimate.most_likely_state,
                left_probabilities=left_estimate.state_probabilities,
                right_probabilities=right_estimate.state_probabilities,
            )
        )
    rows = [
        DiscreteModelComparisonRow(
            model=left_model,
            parameter_count=left.transition_model.parameter_count,
            pseudo_log_likelihood=left.transition_model.pseudo_log_likelihood,
            aic=left.transition_model.aic,
            transition_count=left.transition_summary.transition_count,
        ),
        DiscreteModelComparisonRow(
            model=right_model,
            parameter_count=right.transition_model.parameter_count,
            pseudo_log_likelihood=right.transition_model.pseudo_log_likelihood,
            aic=right.transition_model.aic,
            transition_count=right.transition_summary.transition_count,
        ),
    ]
    better_model = min(rows, key=lambda row: row.aic).model
    sensitive_regions = _build_model_sensitive_regions(differences)
    return DiscreteModelComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=left.taxon_column,
        trait=trait,
        left_model=left_model,
        right_model=right_model,
        better_model=better_model,
        rows=rows,
        node_differences=differences,
        sensitive_region_count=len(sensitive_regions),
        sensitive_regions=sensitive_regions,
    )


def simulate_discrete_stochastic_maps(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    replicates: int = 100,
    seed: int = 0,
) -> StochasticMapCollectionReport:
    """Generate stochastic transition maps from a fitted discrete-state CTMC."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    resolved_model = _resolve_discrete_model_name(model)
    state_order = _resolve_state_order(
        dataset.observed_states,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    fit_report = fit_discrete_mk_model_from_dataset(
        dataset,
        model=resolved_model,
        state_ordering=state_ordering,
        ordered_states=(state_order if state_ordering == "ordered" else None),
    )
    rate_matrix = _rate_matrix_from_transition_rate_rows(
        state_order=fit_report.state_order,
        transition_rate_rows=fit_report.transition_rate_rows,
    )
    root_prior = _resolve_root_prior(
        fit_report.state_order,
        state_counts=dataset.state_counts,
        mode="equal",
        fixed_root_state=None,
        default_root_prior=None,
    )
    return _simulate_stochastic_maps_from_components(
        dataset=dataset,
        model=resolved_model,
        state_order=fit_report.state_order,
        state_ordering=state_ordering,
        ordered_states=(
            fit_report.state_order if state_ordering == "ordered" else []
        ),
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        replicates=replicates,
        seed=seed,
        fit_audit=_stochastic_map_fit_audit(fit_report),
    )


def simulate_discrete_stochastic_maps_from_fit_report(
    fit_report: DiscreteMkFitReport,
    *,
    replicates: int = 100,
    seed: int = 0,
) -> StochasticMapCollectionReport:
    """Generate stochastic maps from one previously fitted discrete Mk surface."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    dataset = load_discrete_dataset(
        fit_report.tree_path,
        fit_report.traits_path,
        trait=fit_report.trait,
        taxon_column=fit_report.taxon_column,
    )
    rate_matrix = _rate_matrix_from_transition_rate_rows(
        state_order=fit_report.state_order,
        transition_rate_rows=fit_report.transition_rate_rows,
    )
    root_prior = _resolve_root_prior(
        fit_report.state_order,
        state_counts=dataset.state_counts,
        mode="equal",
        fixed_root_state=None,
        default_root_prior=None,
    )
    return _simulate_stochastic_maps_from_components(
        dataset=dataset,
        model=fit_report.model,
        state_order=fit_report.state_order,
        state_ordering=fit_report.state_ordering,
        ordered_states=(
            fit_report.state_order if fit_report.state_ordering == "ordered" else []
        ),
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        replicates=replicates,
        seed=seed,
        fit_audit=_stochastic_map_fit_audit(fit_report),
    )


def summarize_discrete_stochastic_maps(
    report: StochasticMapCollectionReport,
) -> StochasticMapSummaryReport:
    """Summarize one stochastic-map collection without regenerating maps."""
    return _summarize_stochastic_map_replicates(
        report.maps,
        simulation_failure_count=len(report.failures),
        expected_transitions=report.fit_audit.allowed_transitions,
    )


def count_discrete_stochastic_map_transitions(
    report: StochasticMapCollectionReport,
) -> StochasticMapTransitionCountReport:
    """Count directional transitions across one stochastic-map collection."""
    summary = summarize_discrete_stochastic_maps(report)
    transition_order = [row.transition for row in summary.rows]
    matrix_rows = [
        StochasticMapTransitionCountMatrixRow(
            replicate_index=replicate.replicate_index,
            total_transition_count=replicate.total_transition_count,
            transition_counts={
                transition: int(replicate.transition_counts.get(transition, 0))
                for transition in transition_order
            },
        )
        for replicate in report.maps
    ]
    branch_keys = [
        (
            history.branch_index,
            history.parent_node,
            history.child_node,
        )
        for history in report.maps[0].branch_histories
    ]
    branch_transition_values: dict[tuple[int, str, str, str], list[int]] = {
        (*branch_key, transition): []
        for branch_key in branch_keys
        for transition in transition_order
    }
    for replicate in report.maps:
        replicate_branch_counts: dict[tuple[int, str, str], dict[str, int]] = {}
        for history in replicate.branch_histories:
            transition_counts = {transition: 0 for transition in transition_order}
            inferred_transitions = [
                f"{event.source_state}->{event.target_state}"
                for event in history.events
            ]
            if not inferred_transitions and len(history.segments) > 1:
                inferred_transitions = [
                    f"{left.state}->{right.state}"
                    for left, right in zip(
                        history.segments,
                        history.segments[1:],
                        strict=False,
                    )
                    if left.state != right.state
                ]
            for transition in inferred_transitions:
                transition_counts[transition] = transition_counts.get(transition, 0) + 1
            replicate_branch_counts[
                (
                    history.branch_index,
                    history.parent_node,
                    history.child_node,
                )
            ] = transition_counts
        for branch_key in branch_keys:
            transition_counts = replicate_branch_counts.get(
                branch_key,
                {transition: 0 for transition in transition_order},
            )
            for transition in transition_order:
                branch_transition_values[(*branch_key, transition)].append(
                    int(transition_counts.get(transition, 0))
                )
    branch_rows: list[StochasticMapBranchTransitionCountRow] = []
    for (
        branch_index,
        parent_node,
        child_node,
        transition,
    ), values in sorted(
        branch_transition_values.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
            item[0][2],
            item[0][3],
        ),
    ):
        sorted_values = sorted(float(value) for value in values)
        branch_rows.append(
            StochasticMapBranchTransitionCountRow(
                branch_index=branch_index,
                parent_node=parent_node,
                child_node=child_node,
                transition=transition,
                mean_count=float(format(sum(values) / max(len(values), 1), ".15g")),
                lower_95_interval=_quantile(sorted_values, 0.025),
                upper_95_interval=_quantile(sorted_values, 0.975),
                minimum_count=min(values, default=0),
                maximum_count=max(values, default=0),
                presence_fraction=float(
                    format(
                        sum(1 for value in values if value > 0) / max(len(values), 1),
                        ".15g",
                    )
                ),
            )
        )
    return StochasticMapTransitionCountReport(
        replicate_count=summary.replicate_count,
        mean_total_transition_count=summary.mean_total_transition_count,
        lower_95_total_transition_count=summary.lower_95_total_transition_count,
        upper_95_total_transition_count=summary.upper_95_total_transition_count,
        transition_order=transition_order,
        matrix_rows=matrix_rows,
        aggregate_rows=summary.rows,
        branch_rows=branch_rows,
        warnings=list(summary.warnings),
    )


def build_biogeographic_interpretation_report(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    compare_model: str | None = None,
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    coding_map: dict[str, str] | None = None,
) -> BiogeographicInterpretationReport:
    """Separate computed biogeographic outputs from downstream interpretation guidance."""
    readiness = assess_geographic_state_analysis_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    report = estimate_ancestral_geographic_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    comparison = (
        compare_discrete_state_models(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            left_model=model,
            right_model=compare_model,
            allowed_states=allowed_states,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
        if compare_model is not None
        else None
    )
    coding_audit = audit_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        coding_map=coding_map,
    )
    root_estimate = next(
        estimate
        for estimate in report.estimates
        if estimate.node == node_signature(load_tree(tree_path).root)
    )
    computed_results = [
        BiogeographicComputedResult(
            label="root_state", value=root_estimate.most_likely_state
        ),
        BiogeographicComputedResult(
            label="transition_count",
            value=str(report.transition_summary.transition_count),
        ),
        BiogeographicComputedResult(
            label="strongly_supported_transition_count",
            value=str(report.transition_summary.strongly_supported_transition_count),
        ),
        BiogeographicComputedResult(
            label="state_count", value=str(len(report.observed_states))
        ),
        BiogeographicComputedResult(label="model", value=report.model),
    ]
    return BiogeographicInterpretationReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        model=model,
        compare_model=compare_model,
        computed_results=computed_results,
        model_sensitive_regions=[]
        if comparison is None
        else comparison.sensitive_regions,
        coding_audit_summary={
            "row_count": coding_audit.row_count,
            "included_row_count": coding_audit.included_row_count,
            "excluded_row_count": coding_audit.excluded_row_count,
        },
        readiness_blockers=readiness.blockers,
        caveats=_build_discrete_evolution_narrative(
            report, comparison=comparison
        ).caveats,
        interpretation_guidance=[
            "computed ancestral regions summarize model-conditioned state histories, not direct evidence of dispersal mechanism",
            "biological interpretation should be restricted to patterns that remain stable across supported models and coding assumptions",
            "sampling gaps, sparse states, and dominant-state bias should be discussed before turning transitions into historical claims",
        ],
    )


def validate_discrete_transition_reference_examples(
    *,
    tolerance: float = 1e-9,
) -> DiscreteTransitionReferenceValidationReport:
    """Validate deterministic discrete-state transition outputs against built-in small references."""
    cases = (
        {
            "label": "toy-geography-er",
            "model": "equal-rates",
            "tree_newick": "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.1);",
            "trait_rows": [
                ("A", "north"),
                ("B", "north"),
                ("C", "south"),
                ("D", "island"),
            ],
            "expected_parameter_count": 1,
            "expected_transition_count": 2,
            "expected_root_state": "north",
            "expected_pseudo_log_likelihood": -7.28950386047299,
            "expected_rates": {
                ("island", "north"): 0.1,
                ("north", "south"): 0.1,
                ("south", "south"): 0.8,
            },
        },
        {
            "label": "toy-geography-sym",
            "model": "symmetric",
            "tree_newick": "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.1);",
            "trait_rows": [
                ("A", "north"),
                ("B", "north"),
                ("C", "south"),
                ("D", "island"),
            ],
            "expected_parameter_count": 3,
            "expected_transition_count": 2,
            "expected_root_state": "north",
            "expected_pseudo_log_likelihood": -6.5047894874944,
            "expected_rates": {
                ("island", "south"): 0.133333333333333,
                ("north", "south"): 0.133333333333333,
                ("south", "island"): 0.1,
            },
        },
        {
            "label": "toy-geography-ard",
            "model": "all-rates-different",
            "tree_newick": "((A:0.1,B:0.1):0.2,(C:0.2,D:0.2):0.1);",
            "trait_rows": [
                ("A", "north"),
                ("B", "north"),
                ("C", "south"),
                ("D", "island"),
            ],
            "expected_parameter_count": 6,
            "expected_transition_count": 2,
            "expected_root_state": "north",
            "expected_pseudo_log_likelihood": -6.2424252230423,
            "expected_rates": {
                ("island", "north"): 0.1,
                ("north", "south"): 0.133333333333333,
                ("south", "island"): 0.133333333333333,
            },
        },
    )
    observations: list[DiscreteTransitionReferenceObservation] = []
    for case in cases:
        with tempfile.TemporaryDirectory(prefix="bijux-discrete-reference-") as tmp_dir:
            tree_path = Path(tmp_dir) / "reference-tree.nwk"
            traits_path = Path(tmp_dir) / "reference-traits.tsv"
            tree_path.write_text(f"{case['tree_newick']}\n", encoding="utf-8")
            traits_path.write_text(
                "taxon\tregion\n"
                + "".join(f"{taxon}\t{state}\n" for taxon, state in case["trait_rows"]),
                encoding="utf-8",
            )
            report = run_discrete_state_transition_model(
                tree_path,
                traits_path,
                trait="region",
                model=str(case["model"]),
            )
        rate_lookup = {
            (row.source_state, target): value
            for row in report.transition_model.transition_matrix
            for target, value in row.target_rates.items()
        }
        rate_rows = [
            DiscreteTransitionReferenceRate(
                source_state=source_state,
                target_state=target_state,
                expected_rate=expected_rate,
                observed_rate=rate_lookup[(source_state, target_state)],
                absolute_delta=abs(
                    rate_lookup[(source_state, target_state)] - expected_rate
                ),
            )
            for (source_state, target_state), expected_rate in sorted(
                case["expected_rates"].items()
            )
        ]
        max_rate_delta = max((row.absolute_delta for row in rate_rows), default=0.0)
        root_state = next(
            estimate.most_likely_state
            for estimate in report.estimates
            if estimate.node == "A|B|C|D"
        )
        passed = (
            report.transition_model.parameter_count == case["expected_parameter_count"]
            and report.transition_summary.transition_count
            == case["expected_transition_count"]
            and root_state == case["expected_root_state"]
            and abs(
                report.transition_model.pseudo_log_likelihood
                - case["expected_pseudo_log_likelihood"]
            )
            <= tolerance
            and max_rate_delta <= tolerance
        )
        observations.append(
            DiscreteTransitionReferenceObservation(
                label=str(case["label"]),
                model=str(case["model"]),
                expected_parameter_count=int(case["expected_parameter_count"]),
                observed_parameter_count=report.transition_model.parameter_count,
                expected_transition_count=int(case["expected_transition_count"]),
                observed_transition_count=report.transition_summary.transition_count,
                expected_root_state=str(case["expected_root_state"]),
                observed_root_state=root_state,
                expected_pseudo_log_likelihood=float(
                    case["expected_pseudo_log_likelihood"]
                ),
                observed_pseudo_log_likelihood=report.transition_model.pseudo_log_likelihood,
                max_rate_delta=max_rate_delta,
                rate_rows=rate_rows,
                passed=passed,
            )
        )
    return DiscreteTransitionReferenceValidationReport(
        case_count=len(observations),
        all_passed=all(observation.passed for observation in observations),
        tolerance=tolerance,
        observations=observations,
    )


def write_node_state_probability_table(
    path: Path, report: DiscreteStateEvolutionReport
) -> Path:
    """Export one deterministic node-probability table for a discrete-state reconstruction."""
    rows = [
        {
            "node": estimate.node,
            "node_name": estimate.node_name or "",
            "is_tip": str(estimate.is_tip).lower(),
            "descendant_taxa": ",".join(estimate.descendant_taxa),
            "most_likely_state": estimate.most_likely_state,
            "state_probabilities": json.dumps(
                estimate.state_probabilities, sort_keys=True
            ),
            "ambiguous": str(estimate.ambiguous).lower(),
        }
        for estimate in report.estimates
    ]
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "is_tip",
            "descendant_taxa",
            "most_likely_state",
            "state_probabilities",
            "ambiguous",
        ],
        rows=rows,
    )


def write_transition_summary_table(
    path: Path, report: DiscreteStateEvolutionReport
) -> Path:
    """Export one branch-by-branch transition summary table."""
    support_by_branch = {
        (row.parent_node, row.child_node): row
        for row in report.transition_summary.support_rows
    }
    rows = [
        {
            "parent_node": event.parent_node,
            "child_node": event.child_node,
            "source_state": event.source_state,
            "target_state": event.target_state,
            "changed": str(event.changed).lower(),
            "support": support_by_branch[(event.parent_node, event.child_node)].support,
            "strongly_supported": str(
                support_by_branch[
                    (event.parent_node, event.child_node)
                ].strongly_supported
            ).lower(),
        }
        for event in report.transition_summary.events
    ]
    return write_taxon_rows(
        path,
        columns=[
            "parent_node",
            "child_node",
            "source_state",
            "target_state",
            "changed",
            "support",
            "strongly_supported",
        ],
        rows=rows,
    )


def write_discrete_model_comparison_table(
    path: Path, report: DiscreteModelComparisonReport
) -> Path:
    """Export one node-wise comparison table across two discrete-state models."""
    rows = [
        {
            "node": difference.node,
            "descendant_taxa": ",".join(difference.descendant_taxa),
            "left_state": difference.left_state,
            "right_state": difference.right_state,
            "differs": str(difference.differs).lower(),
            "left_probabilities": json.dumps(
                difference.left_probabilities, sort_keys=True
            ),
            "right_probabilities": json.dumps(
                difference.right_probabilities, sort_keys=True
            ),
        }
        for difference in report.node_differences
    ]
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "left_state",
            "right_state",
            "differs",
            "left_probabilities",
            "right_probabilities",
        ],
        rows=rows,
    )


def write_stochastic_map_summary_table(
    path: Path, report: StochasticMapSummaryReport
) -> Path:
    """Export one transition-by-transition stochastic-map uncertainty table."""
    rows = [
        {
            "transition": row.transition,
            "mean_count": row.mean_count,
            "lower_95_interval": row.lower_95_interval,
            "upper_95_interval": row.upper_95_interval,
            "minimum_count": row.minimum_count,
            "maximum_count": row.maximum_count,
            "presence_fraction": row.presence_fraction,
        }
        for row in report.rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "transition",
            "mean_count",
            "lower_95_interval",
            "upper_95_interval",
            "minimum_count",
            "maximum_count",
            "presence_fraction",
        ],
        rows=rows,
    )


def write_stochastic_map_transition_count_matrix(
    path: Path, report: StochasticMapTransitionCountReport
) -> Path:
    """Export one countSimmap-style transition matrix with one row per replicate."""
    columns = ["replicate_index", "total_transition_count", *report.transition_order]
    rows = [
        {
            "replicate_index": row.replicate_index,
            "total_transition_count": row.total_transition_count,
            **{
                transition: row.transition_counts.get(transition, 0)
                for transition in report.transition_order
            },
        }
        for row in report.matrix_rows
    ]
    return write_taxon_rows(path, columns=columns, rows=rows)


def write_stochastic_map_aggregate_transition_matrix(
    path: Path, report: StochasticMapTransitionCountReport
) -> Path:
    """Export one mean transition matrix aggregated over a stochastic-map collection."""
    source_states = sorted(
        {
            transition.split("->", 1)[0]
            for transition in report.transition_order
            if "->" in transition
        }
    )
    target_states = sorted(
        {
            transition.split("->", 1)[1]
            for transition in report.transition_order
            if "->" in transition
        }
    )
    mean_lookup = {
        row.transition: row.mean_count
        for row in report.aggregate_rows
    }
    rows = [
        {
            "source_state": source_state,
            **{
                target_state: mean_lookup.get(f"{source_state}->{target_state}", 0.0)
                for target_state in target_states
            },
        }
        for source_state in source_states
    ]
    return write_taxon_rows(
        path,
        columns=["source_state", *target_states],
        rows=rows,
    )


def write_stochastic_map_branch_transition_count_table(
    path: Path, report: StochasticMapTransitionCountReport
) -> Path:
    """Export one per-branch transition-count summary table for a stochastic-map collection."""
    rows = [
        {
            "branch_index": row.branch_index,
            "parent_node": row.parent_node,
            "child_node": row.child_node,
            "transition": row.transition,
            "mean_count": row.mean_count,
            "lower_95_interval": row.lower_95_interval,
            "upper_95_interval": row.upper_95_interval,
            "minimum_count": row.minimum_count,
            "maximum_count": row.maximum_count,
            "presence_fraction": row.presence_fraction,
        }
        for row in report.branch_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "branch_index",
            "parent_node",
            "child_node",
            "transition",
            "mean_count",
            "lower_95_interval",
            "upper_95_interval",
            "minimum_count",
            "maximum_count",
            "presence_fraction",
        ],
        rows=rows,
    )


def write_stochastic_map_state_time_table(
    path: Path, report: StochasticMapSummaryReport
) -> Path:
    """Export one per-state time-in-state summary table for a stochastic-map collection."""
    rows = [
        {
            "state": row.state,
            "mean_time": row.mean_time,
            "lower_95_interval": row.lower_95_interval,
            "upper_95_interval": row.upper_95_interval,
            "minimum_time": row.minimum_time,
            "maximum_time": row.maximum_time,
        }
        for row in report.state_time_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "state",
            "mean_time",
            "lower_95_interval",
            "upper_95_interval",
            "minimum_time",
            "maximum_time",
        ],
        rows=rows,
    )


def write_stochastic_map_branch_occupancy_table(
    path: Path, report: StochasticMapSummaryReport
) -> Path:
    """Export one per-branch state-occupancy summary table for a stochastic-map collection."""
    rows = [
        {
            "branch_index": row.branch_index,
            "parent_node": row.parent_node,
            "child_node": row.child_node,
            "state": row.state,
            "branch_length": row.branch_length,
            "mean_time": row.mean_time,
            "lower_95_interval": row.lower_95_interval,
            "upper_95_interval": row.upper_95_interval,
            "minimum_time": row.minimum_time,
            "maximum_time": row.maximum_time,
            "mean_fraction": row.mean_fraction,
            "presence_fraction": row.presence_fraction,
        }
        for row in report.branch_occupancy_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "branch_index",
            "parent_node",
            "child_node",
            "state",
            "branch_length",
            "mean_time",
            "lower_95_interval",
            "upper_95_interval",
            "minimum_time",
            "maximum_time",
            "mean_fraction",
            "presence_fraction",
        ],
        rows=rows,
    )


def write_stochastic_map_segment_table(
    path: Path, report: StochasticMapCollectionReport
) -> Path:
    """Export one flat branch-state segment table for a stochastic-map collection."""
    rows = [
        {
            "replicate_index": replicate.replicate_index,
            "branch_index": segment.branch_index,
            "parent_node": segment.parent_node,
            "child_node": segment.child_node,
            "state": segment.state,
            "start_time_fraction": segment.start_time_fraction,
            "end_time_fraction": segment.end_time_fraction,
            "duration": segment.duration,
        }
        for replicate in report.maps
        for history in replicate.branch_histories
        for segment in history.segments
    ]
    return write_taxon_rows(
        path,
        columns=[
            "replicate_index",
            "branch_index",
            "parent_node",
            "child_node",
            "state",
            "start_time_fraction",
            "end_time_fraction",
            "duration",
        ],
        rows=rows,
    )


def write_stochastic_map_event_table(
    path: Path, report: StochasticMapCollectionReport
) -> Path:
    """Export one flat transition-event table for a stochastic-map collection."""
    rows = [
        {
            "replicate_index": replicate.replicate_index,
            "branch_index": history.branch_index,
            "parent_node": history.parent_node,
            "child_node": history.child_node,
            "event_index": event_index,
            "source_state": event.source_state,
            "target_state": event.target_state,
            "branch_length": history.branch_length,
            "event_time_fraction": event.event_time_fraction,
            "event_time": float(
                format(history.branch_length * event.event_time_fraction, ".15g")
            ),
        }
        for replicate in report.maps
        for history in replicate.branch_histories
        for event_index, event in enumerate(history.events)
    ]
    return write_taxon_rows(
        path,
        columns=[
            "replicate_index",
            "branch_index",
            "parent_node",
            "child_node",
            "event_index",
            "source_state",
            "target_state",
            "branch_length",
            "event_time_fraction",
            "event_time",
        ],
        rows=rows,
    )


def write_stochastic_map_collection(
    path: Path, report: StochasticMapCollectionReport
) -> Path:
    """Write one stochastic-map collection as JSON."""
    path.write_text(
        json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def load_stochastic_map_collection(path: Path) -> StochasticMapCollectionReport:
    """Load one stochastic-map collection from JSON."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    maps = [
        StochasticMapReplicate(
            replicate_index=replicate["replicate_index"],
            root_state=replicate["root_state"],
            total_transition_count=replicate["total_transition_count"],
            transition_counts=replicate["transition_counts"],
            branch_histories=[
                StochasticMapBranchHistory(
                    branch_index=history["branch_index"],
                    parent_node=history["parent_node"],
                    child_node=history["child_node"],
                    branch_length=history["branch_length"],
                    start_state=history["start_state"],
                    end_state=history["end_state"],
                    event_count=history["event_count"],
                    events=[
                        StochasticMapTransitionEvent(
                            branch_index=event["branch_index"],
                            parent_node=event["parent_node"],
                            child_node=event["child_node"],
                            source_state=event["source_state"],
                            target_state=event["target_state"],
                            event_time_fraction=event["event_time_fraction"],
                        )
                        for event in history["events"]
                    ],
                    segments=[
                        StochasticMapStateSegment(
                            branch_index=segment["branch_index"],
                            parent_node=segment["parent_node"],
                            child_node=segment["child_node"],
                            state=segment["state"],
                            start_time_fraction=segment["start_time_fraction"],
                            end_time_fraction=segment["end_time_fraction"],
                            duration=segment["duration"],
                        )
                        for segment in history.get("segments", [])
                    ],
                )
                for history in replicate["branch_histories"]
            ],
            state_time_totals=replicate.get("state_time_totals", {}),
        )
        for replicate in payload["maps"]
    ]
    summary = StochasticMapSummaryReport(
        replicate_count=payload["summary"]["replicate_count"],
        mean_total_transition_count=payload["summary"]["mean_total_transition_count"],
        lower_95_total_transition_count=payload["summary"][
            "lower_95_total_transition_count"
        ],
        upper_95_total_transition_count=payload["summary"][
            "upper_95_total_transition_count"
        ],
        rows=[
            StochasticMapSummaryRow(
                transition=row["transition"],
                mean_count=row["mean_count"],
                lower_95_interval=row["lower_95_interval"],
                upper_95_interval=row["upper_95_interval"],
                minimum_count=row["minimum_count"],
                maximum_count=row["maximum_count"],
                presence_fraction=row["presence_fraction"],
            )
            for row in payload["summary"]["rows"]
        ],
        state_time_rows=[
            StochasticMapStateTimeRow(
                state=row["state"],
                mean_time=row["mean_time"],
                lower_95_interval=row["lower_95_interval"],
                upper_95_interval=row["upper_95_interval"],
                minimum_time=row["minimum_time"],
                maximum_time=row["maximum_time"],
            )
            for row in payload["summary"].get("state_time_rows", [])
        ],
        branch_occupancy_rows=[
            StochasticMapBranchOccupancyRow(
                branch_index=row["branch_index"],
                parent_node=row["parent_node"],
                child_node=row["child_node"],
                state=row["state"],
                branch_length=row["branch_length"],
                mean_time=row["mean_time"],
                lower_95_interval=row["lower_95_interval"],
                upper_95_interval=row["upper_95_interval"],
                minimum_time=row["minimum_time"],
                maximum_time=row["maximum_time"],
                mean_fraction=row.get("mean_fraction", 0.0),
                presence_fraction=row.get("presence_fraction", 1.0),
            )
            for row in payload["summary"].get("branch_occupancy_rows", [])
        ],
        simulation_failure_count=payload["summary"].get("simulation_failure_count", 0),
        warnings=payload["summary"]["warnings"],
    )
    return StochasticMapCollectionReport(
        tree_path=Path(payload["tree_path"]),
        traits_path=Path(payload["traits_path"]),
        taxon_column=payload["taxon_column"],
        trait=payload["trait"],
        model=payload["model"],
        state_ordering=payload["state_ordering"],
        ordered_states=payload["ordered_states"],
        replicates=payload["replicates"],
        seed=payload["seed"],
        conditioned_on_node_estimates=payload.get(
            "conditioned_on_node_estimates", False
        ),
        fit_audit=StochasticMapModelFitAudit(
            state_order=payload.get("fit_audit", {}).get("state_order", []),
            allowed_transitions=payload.get("fit_audit", {}).get(
                "allowed_transitions", []
            ),
            parameter_count=payload.get("fit_audit", {}).get("parameter_count", 0),
            log_likelihood=payload.get("fit_audit", {}).get("log_likelihood", 0.0),
            aic=payload.get("fit_audit", {}).get("aic", 0.0),
            aicc=payload.get("fit_audit", {}).get("aicc", 0.0),
            overparameterized=payload.get("fit_audit", {}).get(
                "overparameterized", False
            ),
            optimizer_converged=payload.get("fit_audit", {}).get(
                "optimizer_converged", True
            ),
            optimizer_iteration_count=payload.get("fit_audit", {}).get(
                "optimizer_iteration_count", 0
            ),
            optimizer_function_evaluation_count=payload.get("fit_audit", {}).get(
                "optimizer_function_evaluation_count", 0
            ),
            optimizer_hit_lower_parameter_bound=payload.get("fit_audit", {}).get(
                "optimizer_hit_lower_parameter_bound", False
            ),
            optimizer_hit_upper_parameter_bound=payload.get("fit_audit", {}).get(
                "optimizer_hit_upper_parameter_bound", False
            ),
            baseline_model=payload.get("fit_audit", {}).get("baseline_model"),
            baseline_aic=payload.get("fit_audit", {}).get("baseline_aic"),
            baseline_delta_aic=payload.get("fit_audit", {}).get("baseline_delta_aic"),
            preferred_model_by_aic=payload.get("fit_audit", {}).get(
                "preferred_model_by_aic"
            ),
            warnings=payload.get("fit_audit", {}).get("warnings", []),
        ),
        warnings=payload.get("warnings", []),
        maps=maps,
        failures=[
            StochasticMapSimulationFailure(
                replicate_index=row["replicate_index"],
                branch_index=row["branch_index"],
                parent_node=row["parent_node"],
                child_node=row["child_node"],
                source_state=row["source_state"],
                target_state=row["target_state"],
                branch_length=row["branch_length"],
                attempt_count=row["attempt_count"],
                reason=row["reason"],
            )
            for row in payload.get("failures", [])
        ],
        summary=summary,
    )


def render_tree_with_geographic_states(
    tree_path: Path,
    report: DiscreteStateEvolutionReport,
    *,
    out_path: Path,
    layout: str = "phylogram",
) -> TreeRenderResult:
    """Render tip and internal discrete states onto a tree SVG."""
    palette = dict(
        zip(sorted(report.observed_states), _DEFAULT_STATE_COLORS, strict=False)
    )
    internal_annotations = {
        estimate.node: estimate.most_likely_state
        for estimate in report.estimates
        if not estimate.is_tip
    }
    internal_annotation_colors = {
        estimate.node: palette.get(estimate.most_likely_state, "#6d28d9")
        for estimate in report.estimates
        if not estimate.is_tip
    }
    categorical_traits = {
        estimate.node_name: estimate.most_likely_state
        for estimate in report.estimates
        if estimate.is_tip and estimate.node_name is not None
    }
    return render_tree_svg(
        tree_path,
        out_path=out_path,
        layout=layout,
        categorical_traits=categorical_traits,
        internal_annotations=internal_annotations,
        internal_annotation_colors=internal_annotation_colors,
    )


def render_discrete_state_evolution_report(
    *,
    tree_path: Path,
    traits_path: Path,
    trait: str,
    out_path: Path,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    allowed_states: list[str] | None = None,
    compare_model: str | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteEvolutionReportBuildResult:
    """Build a deterministic HTML report for one discrete-state evolution analysis."""
    report = run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    comparison = (
        compare_discrete_state_models(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            left_model=model,
            right_model=compare_model,
            allowed_states=allowed_states,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
        if compare_model is not None
        else None
    )
    render_path = out_path.with_suffix(".svg")
    render_result = render_tree_with_geographic_states(
        tree_path, report, out_path=render_path, layout="phylogram"
    )
    narrative = _build_discrete_evolution_narrative(report, comparison=comparison)
    interpretation = build_biogeographic_interpretation_report(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        compare_model=compare_model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    sections = [
        (
            "discrete-state-summary",
            json.dumps(asdict(narrative), default=str, indent=2, sort_keys=True),
        ),
        (
            "discrete-state-evolution",
            json.dumps(asdict(report), default=str, indent=2, sort_keys=True),
        ),
        (
            "biogeographic-interpretation",
            json.dumps(asdict(interpretation), default=str, indent=2, sort_keys=True),
        ),
        (
            "discrete-state-render",
            json.dumps(asdict(render_result), default=str, indent=2, sort_keys=True),
        ),
    ]
    if comparison is not None:
        sections.append(
            (
                "discrete-state-comparison",
                json.dumps(asdict(comparison), default=str, indent=2, sort_keys=True),
            )
        )
    title = f"Bijux Discrete-State Evolution Report: {trait}"
    machine_manifest = {
        "report_kind": "discrete-state-evolution",
        "title": title,
        "tree_path": str(tree_path),
        "traits_path": str(traits_path),
        "trait": trait,
        "model": model,
        "likelihood_method": report.likelihood_method,
        "state_ordering": report.state_ordering,
        "ordered_states": report.ordered_states,
        "caveat_count": len(narrative.caveats),
        "interpretation_sections": [
            "computed_results",
            "caveats",
            "interpretation_guidance",
        ],
        "rendered_tree": str(render_path),
        "sections": [name for name, _ in sections],
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=machine_manifest,
    )
    return DiscreteEvolutionReportBuildResult(
        output_path=out_path,
        report_kind="discrete-state-evolution",
        title=title,
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        model=model,
        machine_manifest=machine_manifest,
    )
