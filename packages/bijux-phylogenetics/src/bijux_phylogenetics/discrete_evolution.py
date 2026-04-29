from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import random
import tempfile

from bijux_phylogenetics.ancestral.common import load_discrete_dataset, node_descendant_taxa, node_signature
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
class StochasticMapTransitionEvent:
    branch_index: int
    parent_node: str
    child_node: str
    source_state: str
    target_state: str
    event_time_fraction: float


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


@dataclass(slots=True)
class StochasticMapReplicate:
    replicate_index: int
    root_state: str
    total_transition_count: int
    transition_counts: dict[str, int]
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
class StochasticMapSummaryReport:
    replicate_count: int
    mean_total_transition_count: float
    lower_95_total_transition_count: float
    upper_95_total_transition_count: float
    rows: list[StochasticMapSummaryRow]
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
    maps: list[StochasticMapReplicate]
    summary: StochasticMapSummaryReport


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
        return {
            state: float(format(uniform, ".15g"))
            for state in sorted(scores)
        }
    return {
        state: float(format(scores[state] / total, ".15g"))
        for state in sorted(scores)
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


def _state_support(node, states_by_taxon: dict[str, str], state_order: list[str]) -> dict[str, int]:
    counts = {state: 0 for state in state_order}
    for taxon in node_descendant_taxa(node):
        state = states_by_taxon.get(taxon)
        if state is not None:
            counts[state] += 1
    return counts


def _fitch_candidate_sets(tree, states_by_taxon: dict[str, str]) -> dict[str, list[str]]:
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

    return {
        node_signature(node): sorted(downpass(node))
        for node in tree.iter_nodes()
    }


def _best_supported_state(candidate_states: list[str], support_counts: dict[str, int], priority_weights: dict[str, float]) -> str:
    return max(
        sorted(candidate_states),
        key=lambda state: (support_counts.get(state, 0), priority_weights.get(state, 0.0), state),
    )


def _resolve_er_states(tree, candidate_sets: dict[str, list[str]], states_by_taxon: dict[str, str], state_order: list[str]) -> dict[str, str]:
    support_by_node = {
        node_signature(node): _state_support(node, states_by_taxon, state_order)
        for node in tree.iter_nodes()
    }
    resolved: dict[str, str] = {}
    root_signature = node_signature(tree.root)
    root_candidates = candidate_sets[root_signature]
    resolved[root_signature] = _best_supported_state(root_candidates, support_by_node[root_signature], {state: 1.0 for state in state_order})

    def visit(node, parent_state: str) -> None:
        for child in node.children:
            signature = node_signature(child)
            candidates = candidate_sets[signature]
            if parent_state in candidates:
                chosen = parent_state
            else:
                chosen = _best_supported_state(candidates, support_by_node[signature], {state: 1.0 for state in state_order})
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


def _stationary_frequencies(states_by_taxon: dict[str, str], state_order: list[str]) -> dict[str, float]:
    total = len(states_by_taxon)
    if total == 0:
        return {state: 0.0 for state in state_order}
    return {
        state: float(format(sum(1 for value in states_by_taxon.values() if value == state) / total, ".15g"))
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
                if target != source and _transition_allowed(source, target, state_order, state_ordering=state_ordering)
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
        if not _transition_allowed(event.source_state, event.target_state, state_order, state_ordering=state_ordering):
            continue
        if model == "symmetric":
            counts[event.source_state][event.target_state] += 1.0
            if _transition_allowed(event.target_state, event.source_state, state_order, state_ordering=state_ordering):
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
                                        and _transition_allowed(source, candidate, state_order, state_ordering=state_ordering)
                                    ),
                                    1,
                                )
                                if _transition_allowed(source, target, state_order, state_ordering=state_ordering)
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
                                    if _transition_allowed(source, target, state_order, state_ordering=state_ordering) and off_total > 0.0
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
        state: (support_by_node[root_signature].get(state, 0) + 1.0) * root_prior.get(state, 0.0)
        for state in (state_order if state_ordering == "ordered" else candidate_sets[root_signature])
    }
    probabilities_by_node[root_signature] = _normalize_probabilities(root_scores)
    resolved_states[root_signature] = max(sorted(probabilities_by_node[root_signature]), key=lambda state: probabilities_by_node[root_signature][state])

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
                for child_state in (state_order if state_ordering == "ordered" else candidate_sets[child_signature]):
                    transition_support = sum(
                        parent_probabilities[parent_state] * transition_lookup[parent_state][child_state]
                        for parent_state in parent_probabilities
                    )
                    scores[child_state] = transition_support * (child_support.get(child_state, 0) + 1.0)
                probabilities_by_node[child_signature] = _normalize_probabilities(scores)
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
                ambiguous=sum(1 for probability in state_probabilities.values() if probability > 0.0) > 1,
            )
        )
    return estimates


def _root_prior(model: str, stationary: dict[str, float], candidate_states: list[str]) -> dict[str, float]:
    if model == "equal-rates":
        return {state: 1.0 for state in candidate_states}
    return {state: stationary.get(state, 0.0) + 1e-9 for state in candidate_states}


def _pseudo_log_likelihood(estimates: list[NodeStateEstimate], events: list[TransitionEvent], model: TransitionModelReport) -> float:
    estimate_by_node = {estimate.node: estimate for estimate in estimates}
    transition_lookup = _row_lookup(model.transition_matrix)
    root_state = estimates[0].most_likely_state
    log_likelihood = math.log(max(model.root_state_probabilities.get(root_state, 1e-12), 1e-12))
    for event in events:
        probability = transition_lookup[event.source_state][event.target_state]
        log_likelihood += math.log(max(probability, 1e-12))
        child_estimate = estimate_by_node[event.child_node]
        log_likelihood += math.log(max(child_estimate.state_probabilities.get(event.target_state, 1e-12), 1e-12))
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
                    proportion = count_matrix[row.source_state][target_state] / off_total
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


def _summarize_dominant_state_bias(state_counts: dict[str, int]) -> DominantStateBiasReport:
    if not state_counts:
        return DominantStateBiasReport(
            dominant_states=[],
            dominant_fraction=0.0,
            biased=False,
            message=None,
        )
    total = sum(state_counts.values())
    max_count = max(state_counts.values())
    dominant_states = sorted(state for state, count in state_counts.items() if count == max_count)
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
        blockers.append("tree validation failed and the geographic analysis is not safe to interpret")
    if not tree_validation.rooted:
        blockers.append("geographic ancestral-state analysis requires a rooted tree")
    if not coding.valid:
        blockers.append("discrete geographic states contain unsupported labels or coding patterns")
    if any(warning.code == "single-state-dataset" for warning in imbalance.warnings):
        blockers.append("geographic analysis requires at least two observed states after matching taxa to the tree")
    rare_state_count = sum(1 for count in imbalance.state_counts.values() if count < 2)
    if imbalance.state_counts and rare_state_count == len(imbalance.state_counts):
        blockers.append("one or more geographic states are too sparse to estimate transitions credibly")
    if dominant_state_bias.biased:
        blockers.append("observed geographic states are dominated by one state and the sampling is too biased for credible transition inference")

    warnings.extend(tree_validation.warnings)
    warnings.extend(warning.message for warning in imbalance.warnings)
    if dominant_state_bias.message is not None and dominant_state_bias.message not in warnings:
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
        support = 0.0 if total_score <= 0.0 else scores.get(inferred_key, 0.0) / total_score
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
    index = max(0, min(len(sorted_values) - 1, int(round(fraction * (len(sorted_values) - 1)))))
    return float(format(sorted_values[index], ".15g"))


def _build_discrete_evolution_narrative(
    report: DiscreteStateEvolutionReport,
    *,
    comparison: DiscreteModelComparisonReport | None = None,
) -> DiscreteEvolutionNarrative:
    root_state = next(estimate.most_likely_state for estimate in report.estimates if estimate.node == node_signature(load_tree(report.tree_path).root))
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
        right_probability = difference.right_probabilities.get(difference.right_state, 0.0)
        rows.append(
            ModelSensitiveRegionRow(
                node=difference.node,
                descendant_taxa=difference.descendant_taxa,
                left_state=difference.left_state,
                right_state=difference.right_state,
                sensitivity_score=float(format(abs(left_probability - right_probability), ".15g")),
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


def _summarize_stochastic_map_replicates(
    replicates: list[StochasticMapReplicate],
) -> StochasticMapSummaryReport:
    total_counts = sorted(float(replicate.total_transition_count) for replicate in replicates)
    transition_names = sorted(
        {
            transition
            for replicate in replicates
            for transition in replicate.transition_counts
        }
    )
    rows: list[StochasticMapSummaryRow] = []
    for transition in transition_names:
        values = [replicate.transition_counts.get(transition, 0) for replicate in replicates]
        sorted_values = sorted(float(value) for value in values)
        rows.append(
            StochasticMapSummaryRow(
                transition=transition,
                mean_count=float(format(sum(values) / max(len(values), 1), ".15g")),
                lower_95_interval=_quantile(sorted_values, 0.025),
                upper_95_interval=_quantile(sorted_values, 0.975),
                minimum_count=min(values, default=0),
                maximum_count=max(values, default=0),
                presence_fraction=float(format(sum(1 for value in values if value > 0) / max(len(values), 1), ".15g")),
            )
        )
    return StochasticMapSummaryReport(
        replicate_count=len(replicates),
        mean_total_transition_count=float(format(sum(total_counts) / max(len(total_counts), 1), ".15g")),
        lower_95_total_transition_count=_quantile(total_counts, 0.025),
        upper_95_total_transition_count=_quantile(total_counts, 0.975),
        rows=rows,
        warnings=[
            "stochastic maps are simulated from deterministic node estimates and transition weights rather than sampled from a full posterior process"
        ],
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
    table = load_tsv_summary(traits_path) if taxon_column is None else load_taxon_table(traits_path, taxon_column=taxon_column)
    if trait not in table.columns:
        raise AncestralReconstructionError(f"trait table does not contain column '{trait}'")
    tree_taxa = set(tree.tip_names)
    if state_ordering not in {"unordered", "ordered"}:
        raise ValueError(f"unsupported state ordering: {state_ordering}")
    ordered = list(ordered_states or [])
    if ordered and len(set(ordered)) != len(ordered):
        raise AncestralReconstructionError("ordered state vocabulary contains duplicate labels")
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
        invalid_delimiter = next((token for token in _DEFAULT_ALLOWED_STATE_PATTERN_BLOCKLIST if token in raw_state), None)
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
                    code="unordered-state-vocabulary" if state_ordering == "ordered" and ordered else "unsupported-state-label",
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
    table = load_tsv_summary(traits_path) if taxon_column is None else load_taxon_table(traits_path, taxon_column=taxon_column)
    if trait not in table.columns:
        raise AncestralReconstructionError(f"trait table does not contain column '{trait}'")
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
        invalid_delimiter = next((token for token in _DEFAULT_ALLOWED_STATE_PATTERN_BLOCKLIST if token in normalized_state), None)
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
                    issue_code="unsupported-state-label" if state_ordering != "ordered" else "unordered-state-vocabulary",
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
    table = load_tsv_summary(traits_path) if taxon_column is None else load_taxon_table(traits_path, taxon_column=taxon_column)
    if trait not in table.columns:
        raise AncestralReconstructionError(f"trait table does not contain column '{trait}'")
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
    dominant_fraction = max(state_counts.values()) / max(sum(state_counts.values()), 1) if state_counts else 0.0
    if dominant_fraction >= 0.8 and observed_states:
        dominant_states = [state for state, count in state_counts.items() if count == max(state_counts.values())]
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
        raise AncestralReconstructionError("discrete-state evolution input contains unsupported state labels")
    imbalance = detect_state_imbalance_problems(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    if any(warning.code == "single-state-dataset" for warning in imbalance.warnings):
        raise AncestralReconstructionError("discrete-state evolution requires at least two observed states")
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
    er_resolved = _resolve_er_states(dataset.tree, candidate_sets, dataset.states_by_taxon, state_order)
    er_events = _transition_events(dataset.tree, er_resolved)
    count_matrix = _build_transition_count_matrix(
        state_order,
        er_events,
        model=model,
        state_ordering=state_ordering,
    )
    matrix = _fit_transition_matrix(model, state_order, stationary, er_events, state_ordering=state_ordering)
    root_prior = _root_prior(model, stationary, candidate_sets[node_signature(dataset.tree.root)])
    estimates = _estimate_node_states(
        dataset.tree,
        candidate_sets,
        dataset.states_by_taxon,
        state_order,
        matrix,
        root_prior,
        state_ordering=state_ordering,
    )
    resolved_states = {estimate.node: estimate.most_likely_state for estimate in estimates}
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
    strongly_supported_transition_count = sum(1 for row in support_rows if row.strongly_supported)
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
    transition_model.pseudo_log_likelihood = _pseudo_log_likelihood(estimates, events, transition_model)
    transition_model.aic = float(format(2.0 * transition_model.parameter_count - 2.0 * transition_model.pseudo_log_likelihood, ".15g"))
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
        strongly_supported_transition_counts=dict(sorted(strongly_supported_transition_counts.items())),
        support_rows=support_rows,
        events=events,
    )
    warnings = list(dataset.warnings)
    warnings.extend(warning.message for warning in imbalance.warnings)
    if instability.unstable:
        warnings.append("sparse-state transition estimates may be unstable for one or more source-target paths")
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
            "geographic state analysis is inappropriate: " + "; ".join(readiness.blockers)
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
                differs=left_estimate.most_likely_state != right_estimate.most_likely_state,
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
    """Generate approximate stochastic transition maps from deterministic node-state estimates."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
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
    tree = load_tree(tree_path)
    estimate_by_node = {estimate.node: estimate for estimate in report.estimates}
    transition_lookup = _row_lookup(report.transition_model.transition_matrix)
    branch_rows = [
        (index, event.parent_node, event.child_node, next(node for node in tree.iter_nodes() if node_signature(node) == event.child_node))
        for index, event in enumerate(report.transition_summary.events)
    ]
    randomizer = random.Random(seed)
    maps: list[StochasticMapReplicate] = []
    for replicate_index in range(replicates):
        root_state = _sample_state(
            next(estimate.state_probabilities for estimate in report.estimates if estimate.node == node_signature(tree.root)),
            randomizer,
        )
        node_states = {node_signature(tree.root): root_state}
        branch_histories: list[StochasticMapBranchHistory] = []
        transition_counts: dict[str, int] = {}
        total_transition_count = 0
        for branch_index, parent_node, child_node, child in branch_rows:
            parent_state = node_states[parent_node]
            child_estimate = estimate_by_node[child_node]
            branch_length = float(child.branch_length or 0.0)
            start_state = parent_state
            end_state = _sample_state(child_estimate.state_probabilities, randomizer)
            max_event_count = max(1, int(math.ceil(max(branch_length, 0.0) * 2.0)))
            provisional_event_count = 0 if start_state == end_state else randomizer.randint(1, max_event_count)
            events: list[StochasticMapTransitionEvent] = []
            current_state = start_state
            for _ in range(provisional_event_count):
                next_state = _sample_transition_target(current_state, transition_lookup, randomizer)
                if next_state == current_state:
                    continue
                event = StochasticMapTransitionEvent(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    source_state=current_state,
                    target_state=next_state,
                    event_time_fraction=float(format(randomizer.random(), ".15g")),
                )
                events.append(event)
                transition = f"{current_state}->{next_state}"
                transition_counts[transition] = transition_counts.get(transition, 0) + 1
                total_transition_count += 1
                current_state = next_state
            if current_state != end_state:
                events.append(
                    StochasticMapTransitionEvent(
                        branch_index=branch_index,
                        parent_node=parent_node,
                        child_node=child_node,
                        source_state=current_state,
                        target_state=end_state,
                        event_time_fraction=1.0,
                    )
                )
                transition = f"{current_state}->{end_state}"
                transition_counts[transition] = transition_counts.get(transition, 0) + 1
                total_transition_count += 1
            ordered_events = sorted(events, key=lambda event: event.event_time_fraction)
            node_states[child_node] = end_state
            branch_histories.append(
                StochasticMapBranchHistory(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    branch_length=branch_length,
                    start_state=start_state,
                    end_state=end_state,
                    event_count=len(ordered_events),
                    events=ordered_events,
                )
            )
        maps.append(
            StochasticMapReplicate(
                replicate_index=replicate_index,
                root_state=root_state,
                total_transition_count=total_transition_count,
                transition_counts=dict(sorted(transition_counts.items())),
                branch_histories=branch_histories,
            )
        )
    return StochasticMapCollectionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=report.taxon_column,
        trait=trait,
        model=model,
        state_ordering=report.state_ordering,
        ordered_states=report.ordered_states,
        replicates=replicates,
        seed=seed,
        conditioned_on_node_estimates=True,
        maps=maps,
        summary=_summarize_stochastic_map_replicates(maps),
    )


def summarize_discrete_stochastic_maps(
    report: StochasticMapCollectionReport,
) -> StochasticMapSummaryReport:
    """Summarize one stochastic-map collection without regenerating maps."""
    return _summarize_stochastic_map_replicates(report.maps)


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
            "trait_rows": [("A", "north"), ("B", "north"), ("C", "south"), ("D", "island")],
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
            "trait_rows": [("A", "north"), ("B", "north"), ("C", "south"), ("D", "island")],
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
            "trait_rows": [("A", "north"), ("B", "north"), ("C", "south"), ("D", "island")],
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
                absolute_delta=abs(rate_lookup[(source_state, target_state)] - expected_rate),
            )
            for (source_state, target_state), expected_rate in sorted(case["expected_rates"].items())
        ]
        max_rate_delta = max((row.absolute_delta for row in rate_rows), default=0.0)
        root_state = next(estimate.most_likely_state for estimate in report.estimates if estimate.node == "A|B|C|D")
        passed = (
            report.transition_model.parameter_count == case["expected_parameter_count"]
            and report.transition_summary.transition_count == case["expected_transition_count"]
            and root_state == case["expected_root_state"]
            and abs(report.transition_model.pseudo_log_likelihood - case["expected_pseudo_log_likelihood"]) <= tolerance
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
                expected_pseudo_log_likelihood=float(case["expected_pseudo_log_likelihood"]),
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


def write_node_state_probability_table(path: Path, report: DiscreteStateEvolutionReport) -> Path:
    """Export one deterministic node-probability table for a discrete-state reconstruction."""
    rows = [
        {
            "node": estimate.node,
            "node_name": estimate.node_name or "",
            "is_tip": str(estimate.is_tip).lower(),
            "descendant_taxa": ",".join(estimate.descendant_taxa),
            "most_likely_state": estimate.most_likely_state,
            "state_probabilities": json.dumps(estimate.state_probabilities, sort_keys=True),
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


def write_transition_summary_table(path: Path, report: DiscreteStateEvolutionReport) -> Path:
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
            "strongly_supported": str(support_by_branch[(event.parent_node, event.child_node)].strongly_supported).lower(),
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


def write_discrete_model_comparison_table(path: Path, report: DiscreteModelComparisonReport) -> Path:
    """Export one node-wise comparison table across two discrete-state models."""
    rows = [
        {
            "node": difference.node,
            "descendant_taxa": ",".join(difference.descendant_taxa),
            "left_state": difference.left_state,
            "right_state": difference.right_state,
            "differs": str(difference.differs).lower(),
            "left_probabilities": json.dumps(difference.left_probabilities, sort_keys=True),
            "right_probabilities": json.dumps(difference.right_probabilities, sort_keys=True),
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


def write_stochastic_map_summary_table(path: Path, report: StochasticMapSummaryReport) -> Path:
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


def write_stochastic_map_collection(path: Path, report: StochasticMapCollectionReport) -> Path:
    """Write one stochastic-map collection as JSON."""
    path.write_text(json.dumps(asdict(report), default=str, indent=2, sort_keys=True) + "\n", encoding="utf-8")
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
                )
                for history in replicate["branch_histories"]
            ],
        )
        for replicate in payload["maps"]
    ]
    summary = StochasticMapSummaryReport(
        replicate_count=payload["summary"]["replicate_count"],
        mean_total_transition_count=payload["summary"]["mean_total_transition_count"],
        lower_95_total_transition_count=payload["summary"]["lower_95_total_transition_count"],
        upper_95_total_transition_count=payload["summary"]["upper_95_total_transition_count"],
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
        conditioned_on_node_estimates=payload["conditioned_on_node_estimates"],
        maps=maps,
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
    palette = {
        state: color
        for state, color in zip(sorted(report.observed_states), _DEFAULT_STATE_COLORS, strict=False)
    }
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
    render_result = render_tree_with_geographic_states(tree_path, report, out_path=render_path, layout="phylogram")
    narrative = _build_discrete_evolution_narrative(report, comparison=comparison)
    sections = [
        ("discrete-state-summary", json.dumps(asdict(narrative), default=str, indent=2, sort_keys=True)),
        ("discrete-state-evolution", json.dumps(asdict(report), default=str, indent=2, sort_keys=True)),
        ("discrete-state-render", json.dumps(asdict(render_result), default=str, indent=2, sort_keys=True)),
    ]
    if comparison is not None:
        sections.append(("discrete-state-comparison", json.dumps(asdict(comparison), default=str, indent=2, sort_keys=True)))
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
        "rendered_tree": str(render_path),
        "sections": [name for name, _ in sections],
    }
    write_html_report(title=title, sections=sections, out_path=out_path, embedded_json=machine_manifest)
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
