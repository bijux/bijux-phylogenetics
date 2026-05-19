from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    load_discrete_dataset,
    node_descendant_taxa,
    node_signature,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)
from bijux_phylogenetics.core.metadata import load_taxon_table
from bijux_phylogenetics.core.traits import load_tsv_summary
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree

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
    if model == "meristic":
        _resolve_discrete_model_name(model)
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
