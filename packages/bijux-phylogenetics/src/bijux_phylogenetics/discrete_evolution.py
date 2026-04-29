from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import load_discrete_dataset, node_descendant_taxa, node_signature
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.core.traits import load_tsv_summary
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
    valid: bool
    issues: list[StateCodingIssue]
    observed_states: list[str]
    usable_taxa: list[str]


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
class TransitionModelReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    state_order: list[str]
    parameter_count: int
    pseudo_log_likelihood: float
    aic: float
    stationary_frequencies: dict[str, float]
    transition_matrix: list[TransitionRateRow]
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
class TransitionSummaryReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    branch_count: int
    transition_count: int
    transition_counts: dict[str, int]
    events: list[TransitionEvent]


@dataclass(slots=True)
class DiscreteStateEvolutionReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    analysis_tree_newick: str
    taxon_count: int
    observed_states: list[str]
    state_counts: dict[str, int]
    coding_validation: StateCodingValidationReport
    imbalance: StateImbalanceReport
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


def _fit_transition_matrix(model: str, state_order: list[str], stationary: dict[str, float], er_events: list[TransitionEvent]) -> list[TransitionRateRow]:
    change_mass = 0.2
    stay_mass = 0.8
    if model == "equal-rates":
        off_diagonal = change_mass / max(len(state_order) - 1, 1)
        return [
            TransitionRateRow(
                source_state=source,
                target_rates={
                    target: float(format(stay_mass if source == target else off_diagonal, ".15g"))
                    for target in state_order
                },
            )
            for source in state_order
        ]

    counts = {
        source: {target: 1.0 for target in state_order if target != source}
        for source in state_order
    }
    for event in er_events:
        if event.source_state != event.target_state:
            counts[event.source_state][event.target_state] += 1.0
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
                        else float(format(change_mass * (counts[source][target] / off_total), ".15g"))
                    )
                    for target in state_order
                },
            )
        )
    return rows


def _row_lookup(rows: list[TransitionRateRow]) -> dict[str, dict[str, float]]:
    return {row.source_state: row.target_rates for row in rows}


def _estimate_node_states(tree, candidate_sets: dict[str, list[str]], states_by_taxon: dict[str, str], state_order: list[str], matrix: list[TransitionRateRow], root_prior: dict[str, float]) -> list[NodeStateEstimate]:
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
        for state in candidate_sets[root_signature]
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
                for child_state in candidate_sets[child_signature]:
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


def validate_discrete_state_coding(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    allowed_states: list[str] | None = None,
) -> StateCodingValidationReport:
    """Detect impossible or unsupported discrete-state labels."""
    tree = load_tree(tree_path)
    table = load_tsv_summary(traits_path) if taxon_column is None else load_taxon_table(traits_path, taxon_column=taxon_column)
    if trait not in table.columns:
        raise AncestralReconstructionError(f"trait table does not contain column '{trait}'")
    tree_taxa = set(tree.tip_names)
    allowed = sorted(allowed_states or [])
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
                    code="unsupported-state-label",
                    message="state label is not present in the allowed state vocabulary",
                )
            )
            continue
        observed_states.add(raw_state)
        usable_taxa.append(taxon)
    return StateCodingValidationReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        allowed_states=allowed,
        valid=not issues,
        issues=issues,
        observed_states=sorted(observed_states),
        usable_taxa=sorted(usable_taxa),
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
) -> DiscreteStateEvolutionReport:
    """Run a deterministic discrete-state evolution workflow on one tree and trait."""
    if model not in {"equal-rates", "all-rates-different"}:
        raise ValueError(f"unsupported discrete-state model: {model}")
    coding = validate_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_states,
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
    state_order = sorted(dataset.observed_states)
    candidate_sets = _fitch_candidate_sets(dataset.tree, dataset.states_by_taxon)
    stationary = _stationary_frequencies(dataset.states_by_taxon, state_order)
    er_resolved = _resolve_er_states(dataset.tree, candidate_sets, dataset.states_by_taxon, state_order)
    er_events = _transition_events(dataset.tree, er_resolved)
    matrix = _fit_transition_matrix(model, state_order, stationary, er_events)
    root_prior = _root_prior(model, stationary, candidate_sets[node_signature(dataset.tree.root)])
    estimates = _estimate_node_states(dataset.tree, candidate_sets, dataset.states_by_taxon, state_order, matrix, root_prior)
    resolved_states = {estimate.node: estimate.most_likely_state for estimate in estimates}
    events = _transition_events(dataset.tree, resolved_states)
    transition_counts: dict[str, int] = {}
    for event in events:
        key = f"{event.source_state}->{event.target_state}"
        transition_counts[key] = transition_counts.get(key, 0) + 1
    branch_count = len(events)
    transition_count = sum(1 for event in events if event.changed)
    transition_model = TransitionModelReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=model,
        state_order=state_order,
        parameter_count=1 if model == "equal-rates" else len(state_order) * max(len(state_order) - 1, 0),
        pseudo_log_likelihood=0.0,
        aic=0.0,
        stationary_frequencies=stationary,
        transition_matrix=matrix,
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
        transition_counts=dict(sorted(transition_counts.items())),
        events=events,
    )
    warnings = list(dataset.warnings)
    warnings.extend(warning.message for warning in imbalance.warnings)
    return DiscreteStateEvolutionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=model,
        analysis_tree_newick=dumps_newick(dataset.tree),
        taxon_count=len(dataset.taxa),
        observed_states=state_order,
        state_counts=dataset.state_counts,
        coding_validation=coding,
        imbalance=imbalance,
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
) -> DiscreteStateEvolutionReport:
    """Estimate ancestral geographic states over a rooted tree."""
    return run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_states=allowed_states,
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
) -> DiscreteModelComparisonReport:
    """Compare discrete-state reconstructions across two supported models."""
    left = run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=left_model,
        allowed_states=allowed_states,
    )
    right = run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=right_model,
        allowed_states=allowed_states,
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
) -> DiscreteEvolutionReportBuildResult:
    """Build a deterministic HTML report for one discrete-state evolution analysis."""
    report = run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_states=allowed_states,
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
        )
        if compare_model is not None
        else None
    )
    render_path = out_path.with_suffix(".svg")
    render_result = render_tree_with_geographic_states(tree_path, report, out_path=render_path, layout="phylogram")
    sections = [
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
