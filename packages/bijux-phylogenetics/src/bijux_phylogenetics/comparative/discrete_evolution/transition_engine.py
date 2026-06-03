from __future__ import annotations

import math

from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature

from .models import (
    NodeStateEstimate,
    TransitionEvent,
    TransitionModelReport,
    TransitionRateRow,
    TransitionRateUncertaintyReport,
    TransitionRateUncertaintyRow,
    TransitionSupportRow,
)


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
