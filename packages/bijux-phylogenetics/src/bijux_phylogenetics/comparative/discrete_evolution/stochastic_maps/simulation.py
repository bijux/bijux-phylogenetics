from __future__ import annotations

import random

import numpy

from bijux_phylogenetics.ancestral.common import load_discrete_dataset, node_signature
from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    branch_length as _branch_length,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    transition_probability_matrix as _transition_probability_matrix,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_root_prior as _resolve_root_prior,
)
from bijux_phylogenetics.comparative.discrete_mk import (
    DiscreteMkFitReport,
    fit_discrete_mk_model_from_dataset,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from ..transition_engine import _resolve_state_order
from .models import (
    StochasticMapBranchHistory,
    StochasticMapCollectionReport,
    StochasticMapModelFitAudit,
    StochasticMapReplicate,
    StochasticMapSimulationFailure,
    StochasticMapStateSegment,
    StochasticMapTransitionEvent,
)
from .summary import _stochastic_map_warning_union, _summarize_stochastic_map_replicates


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
            node_states[child_signature] = state_order[
                _sample_index(child_weights, rng)
            ]
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
        state_time_totals = dict.fromkeys(state_order, 0.0)
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
                state_time_totals = dict.fromkeys(state_order, 0.0)
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


def simulate_discrete_stochastic_maps(
    tree_path,
    traits_path,
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
        ordered_states=(fit_report.state_order if state_ordering == "ordered" else []),
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
