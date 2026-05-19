from __future__ import annotations

from dataclasses import asdict, dataclass
from html import escape
import json
import math
from pathlib import Path
import random

import numpy

from bijux_phylogenetics.ancestral.common import load_discrete_dataset, node_signature
from bijux_phylogenetics.ancestral.discrete import _branch_length, _transition_probability_matrix
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
    resolve_root_prior as _resolve_root_prior,
)
from bijux_phylogenetics.comparative.discrete_mk import (
    DiscreteMkFitReport,
    fit_discrete_mk_model_from_dataset,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.svg import render_tree_svg

from .core import (
    _DEFAULT_STATE_COLORS,
    _quantile,
    _resolve_state_order,
    AncestralReconstructionError,
)


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
        {state for replicate in replicates for state in replicate.state_time_totals}
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
            state_times = dict.fromkeys(state_names, 0.0)
            for segment in history.segments:
                state_times[segment.state] = state_times.get(
                    segment.state, 0.0
                ) + float(segment.duration)
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
    state_order = list(report.fit_audit.state_order)
    if not state_order:
        state_order = sorted(
            {
                state
                for replicate in report.maps
                for state in replicate.state_time_totals
            }
        )
    transition_order = [
        f"{source_state}->{target_state}"
        for source_state in state_order
        for target_state in state_order
    ]
    aggregate_lookup = {row.transition: row for row in summary.rows}
    aggregate_rows = [
        aggregate_lookup.get(
            transition,
            StochasticMapSummaryRow(
                transition=transition,
                mean_count=0.0,
                lower_95_interval=0.0,
                upper_95_interval=0.0,
                minimum_count=0,
                maximum_count=0,
                presence_fraction=0.0,
            ),
        )
        for transition in transition_order
    ]
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
            transition_counts = dict.fromkeys(transition_order, 0)
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
                dict.fromkeys(transition_order, 0),
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
        aggregate_rows=aggregate_rows,
        branch_rows=branch_rows,
        warnings=list(summary.warnings),
    )


def _binary_entropy(probability: float) -> float:
    if probability <= 0.0 or probability >= 1.0:
        return 0.0
    return float(
        format(
            -(probability * math.log2(probability))
            - ((1.0 - probability) * math.log2(1.0 - probability)),
            ".15g",
        )
    )


def _normalize_probability_interval(value: float, branch_length: float) -> float:
    if branch_length <= 0.0:
        return 0.0
    return float(format(value / branch_length, ".15g"))


def _hex_to_rgb(color: str) -> tuple[int, int, int]:
    token = color.removeprefix("#")
    return (
        int(token[0:2], 16),
        int(token[2:4], 16),
        int(token[4:6], 16),
    )


def _blend_density_color(
    probability: float,
    *,
    low_color: str = _DEFAULT_STATE_COLORS[0],
    high_color: str = _DEFAULT_STATE_COLORS[1],
) -> str:
    bounded = min(max(probability, 0.0), 1.0)
    low_red, low_green, low_blue = _hex_to_rgb(low_color)
    high_red, high_green, high_blue = _hex_to_rgb(high_color)
    red = round(low_red + (high_red - low_red) * bounded)
    green = round(low_green + (high_green - low_green) * bounded)
    blue = round(low_blue + (high_blue - low_blue) * bounded)
    return f"#{red:02x}{green:02x}{blue:02x}"


def _tree_branch_geometry(
    report: StochasticMapCollectionReport,
) -> list[tuple[int, str, str, float, float, float]]:
    tree = load_tree(report.tree_path)
    root_candidates = {
        history.parent_node for history in report.maps[0].branch_histories
    } - {history.child_node for history in report.maps[0].branch_histories}
    if len(root_candidates) == 1:
        kept_taxa = next(iter(root_candidates)).split("|")
        if sorted(kept_taxa) != sorted(tree.tip_names):
            tree, _ = prune_tree_to_requested_taxa(report.tree_path, kept_taxa)
    expected_branch_indexes = {
        (
            history.parent_node,
            history.child_node,
        ): history.branch_index
        for history in report.maps[0].branch_histories
    }
    rows: list[tuple[int, str, str, float, float, float]] = []

    def visit(node, parent_depth: float) -> None:
        for child in node.children:
            branch_length = float(child.branch_length or 0.0)
            parent_node = node_signature(node)
            child_node = node_signature(child)
            branch_index = expected_branch_indexes.get((parent_node, child_node))
            if branch_index is not None:
                rows.append(
                    (
                        branch_index,
                        parent_node,
                        child_node,
                        branch_length,
                        float(format(parent_depth, ".15g")),
                        float(format(parent_depth + branch_length, ".15g")),
                    )
                )
            visit(child, parent_depth + branch_length)

    visit(tree.root, 0.0)
    return rows


def summarize_discrete_stochastic_map_density(
    report: StochasticMapCollectionReport,
    *,
    resolution: int = 100,
    focal_state: str | None = None,
) -> StochasticMapDensityReport:
    """Summarize one stochastic-map collection as branch probability density."""
    if resolution < 1:
        raise ValueError(f"resolution must be at least 1, got {resolution}")
    summary = summarize_discrete_stochastic_maps(report)
    branch_state_rows = [
        StochasticMapBranchProbabilityRow(
            branch_index=row.branch_index,
            parent_node=row.parent_node,
            child_node=row.child_node,
            state=row.state,
            branch_length=row.branch_length,
            mean_probability=_normalize_probability_interval(
                row.mean_time, row.branch_length
            ),
            lower_95_probability=_normalize_probability_interval(
                row.lower_95_interval, row.branch_length
            ),
            upper_95_probability=_normalize_probability_interval(
                row.upper_95_interval, row.branch_length
            ),
            minimum_probability=_normalize_probability_interval(
                row.minimum_time, row.branch_length
            ),
            maximum_probability=_normalize_probability_interval(
                row.maximum_time, row.branch_length
            ),
            presence_fraction=row.presence_fraction,
        )
        for row in summary.branch_occupancy_rows
    ]
    observed_state_order = list(dict.fromkeys(row.state for row in branch_state_rows))
    declared_state_order = list(report.fit_audit.state_order)
    if declared_state_order:
        state_order = [
            state for state in declared_state_order if state in observed_state_order
        ]
    else:
        state_order = []
    if not state_order:
        state_order = observed_state_order
    warnings = list(summary.warnings)
    resolved_focal_state = focal_state
    baseline_state: str | None = None
    if resolved_focal_state is None:
        if len(state_order) == 2:
            baseline_state = state_order[0]
            resolved_focal_state = state_order[1]
        elif len(state_order) > 2:
            warnings.append(
                "branch-state probability summaries are available for multistate collections, but density slices require one explicit focal state"
            )
    else:
        allowed_focal_states = declared_state_order or state_order
        if resolved_focal_state not in allowed_focal_states:
            raise ValueError(
                f"focal_state '{resolved_focal_state}' is not present in the stochastic-map state order"
            )
        if len(state_order) == 2:
            baseline_state = next(
                state for state in state_order if state != resolved_focal_state
            )
    branch_rows: list[StochasticMapDensityBranchRow] = []
    density_rows: list[StochasticMapDensitySliceRow] = []
    branch_geometry = _tree_branch_geometry(report)
    total_tree_depth = max(
        (end_depth for _, _, _, _, _, end_depth in branch_geometry),
        default=0.0,
    )
    if resolved_focal_state is None:
        return StochasticMapDensityReport(
            replicate_count=summary.replicate_count,
            resolution=resolution,
            total_tree_depth=total_tree_depth,
            state_order=state_order,
            focal_state=None,
            baseline_state=baseline_state,
            branch_state_rows=branch_state_rows,
            density_rows=density_rows,
            branch_rows=branch_rows,
            warnings=warnings,
        )
    branch_histories_by_replicate = [
        {history.branch_index: history for history in replicate.branch_histories}
        for replicate in report.maps
    ]
    steps = [
        (float(step) / float(resolution)) * total_tree_depth
        for step in range(resolution + 1)
    ]
    for (
        branch_index,
        parent_node,
        child_node,
        branch_length,
        start_depth,
        end_depth,
    ) in branch_geometry:
        boundaries = [start_depth]
        boundaries.extend(step for step in steps if start_depth < step < end_depth)
        boundaries.append(end_depth)
        branch_slice_rows: list[StochasticMapDensitySliceRow] = []
        if branch_length <= 0.0:
            branch_rows.append(
                StochasticMapDensityBranchRow(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    branch_length=branch_length,
                    focal_state=resolved_focal_state,
                    baseline_state=baseline_state,
                    mean_posterior_probability=0.0,
                    minimum_posterior_probability=0.0,
                    maximum_posterior_probability=0.0,
                    uncertainty=0.0,
                    slice_count=0,
                )
            )
            continue
        for slice_index, (slice_start, slice_end) in enumerate(
            zip(boundaries, boundaries[1:], strict=False)
        ):
            slice_length = slice_end - slice_start
            if slice_length <= 0.0:
                continue
            slice_start_local = slice_start - start_depth
            slice_end_local = slice_end - start_depth
            replicate_probabilities: list[float] = []
            for branch_histories in branch_histories_by_replicate:
                history = branch_histories[branch_index]
                focal_duration = 0.0
                for segment in history.segments:
                    segment_start = segment.start_time_fraction * branch_length
                    segment_end = segment.end_time_fraction * branch_length
                    overlap = min(segment_end, slice_end_local) - max(
                        segment_start, slice_start_local
                    )
                    if overlap <= 0.0 or segment.state != resolved_focal_state:
                        continue
                    focal_duration += overlap
                replicate_probabilities.append(focal_duration / slice_length)
            posterior_probability = float(
                format(
                    sum(replicate_probabilities) / max(len(replicate_probabilities), 1),
                    ".15g",
                )
            )
            branch_slice_rows.append(
                StochasticMapDensitySliceRow(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    branch_length=branch_length,
                    slice_index=slice_index,
                    start_depth=float(format(slice_start, ".15g")),
                    end_depth=float(format(slice_end, ".15g")),
                    start_time_fraction=float(
                        format(slice_start_local / branch_length, ".15g")
                    ),
                    end_time_fraction=float(
                        format(slice_end_local / branch_length, ".15g")
                    ),
                    posterior_probability=posterior_probability,
                    posterior_uncertainty=_binary_entropy(posterior_probability),
                )
            )
        density_rows.extend(branch_slice_rows)
        if not branch_slice_rows:
            branch_rows.append(
                StochasticMapDensityBranchRow(
                    branch_index=branch_index,
                    parent_node=parent_node,
                    child_node=child_node,
                    branch_length=branch_length,
                    focal_state=resolved_focal_state,
                    baseline_state=baseline_state,
                    mean_posterior_probability=0.0,
                    minimum_posterior_probability=0.0,
                    maximum_posterior_probability=0.0,
                    uncertainty=0.0,
                    slice_count=0,
                )
            )
            continue
        weighted_probability_total = 0.0
        weighted_uncertainty_total = 0.0
        probability_values: list[float] = []
        for row in branch_slice_rows:
            slice_length = row.end_depth - row.start_depth
            weighted_probability_total += row.posterior_probability * slice_length
            weighted_uncertainty_total += row.posterior_uncertainty * slice_length
            probability_values.append(row.posterior_probability)
        branch_rows.append(
            StochasticMapDensityBranchRow(
                branch_index=branch_index,
                parent_node=parent_node,
                child_node=child_node,
                branch_length=branch_length,
                focal_state=resolved_focal_state,
                baseline_state=baseline_state,
                mean_posterior_probability=float(
                    format(weighted_probability_total / branch_length, ".15g")
                ),
                minimum_posterior_probability=min(probability_values, default=0.0),
                maximum_posterior_probability=max(probability_values, default=0.0),
                uncertainty=float(
                    format(weighted_uncertainty_total / branch_length, ".15g")
                ),
                slice_count=len(branch_slice_rows),
            )
        )
    return StochasticMapDensityReport(
        replicate_count=summary.replicate_count,
        resolution=resolution,
        total_tree_depth=total_tree_depth,
        state_order=state_order,
        focal_state=resolved_focal_state,
        baseline_state=baseline_state,
        branch_state_rows=branch_state_rows,
        density_rows=density_rows,
        branch_rows=branch_rows,
        warnings=warnings,
    )


def _write_stochastic_map_density_html(
    *,
    report: StochasticMapDensityReport,
    svg_path: Path,
    out_path: Path,
    layout: str,
) -> None:
    svg_markup = svg_path.read_text(encoding="utf-8")
    low_color = _DEFAULT_STATE_COLORS[0]
    high_color = _DEFAULT_STATE_COLORS[1]
    summary_cards = [
        ("focal state", report.focal_state or ""),
        ("baseline state", report.baseline_state or "complement"),
        ("resolution", str(report.resolution)),
        ("replicates", str(report.replicate_count)),
        ("branch density rows", str(len(report.branch_rows))),
    ]
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(f"Bijux Stochastic Density Map: {report.focal_state or 'state density'}")}</title>
  <style>
    :root {{
      color-scheme: light;
      --ink: #1b1f24;
      --bg: #f8fafc;
      --panel: #ffffff;
      --rule: #d6dee8;
      --accent: #0f766e;
      --mono: "SFMono-Regular", "SF Mono", Consolas, monospace;
    }}
    body {{
      margin: 0;
      padding: 2rem;
      background: linear-gradient(180deg, #eef6f4 0%, var(--bg) 100%);
      color: var(--ink);
      font: 16px/1.5 "Iowan Old Style", "Palatino Linotype", serif;
    }}
    main {{
      max-width: 1200px;
      margin: 0 auto;
      display: grid;
      gap: 1.5rem;
    }}
    section {{
      background: var(--panel);
      border: 1px solid var(--rule);
      border-radius: 18px;
      padding: 1.5rem;
      box-shadow: 0 20px 60px rgba(15, 118, 110, 0.08);
    }}
    h1, h2 {{
      font-family: "Avenir Next", "Segoe UI", sans-serif;
    }}
    h1 {{
      margin: 0;
      color: var(--accent);
    }}
    .summary-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 0.75rem;
    }}
    .summary-card {{
      background: #f6fffd;
      border: 1px solid rgba(15, 118, 110, 0.18);
      border-radius: 14px;
      padding: 0.9rem 1rem;
    }}
    .summary-card dt {{
      margin: 0;
      font: 600 0.82rem/1.2 "Avenir Next", "Segoe UI", sans-serif;
      letter-spacing: 0.03em;
      text-transform: uppercase;
      color: #476b67;
    }}
    .summary-card dd {{
      margin: 0.35rem 0 0;
      font: 700 1.1rem/1.2 var(--mono);
    }}
    .figure-shell svg {{
      width: 100%;
      height: auto;
      display: block;
    }}
    .legend-bar {{
      width: min(420px, 100%);
      height: 18px;
      border-radius: 999px;
      border: 1px solid rgba(15, 23, 42, 0.14);
      background: linear-gradient(90deg, {low_color} 0%, {high_color} 100%);
    }}
    .legend-labels {{
      display: flex;
      justify-content: space-between;
      font-family: "Avenir Next", "Segoe UI", sans-serif;
      font-size: 0.9rem;
      margin-top: 0.35rem;
    }}
  </style>
</head>
<body>
  <main>
    <section>
      <h1>{escape(report.focal_state or "stochastic density map")}</h1>
      <p>Reviewer-facing stochastic-map density artifact with branch colors scaled by mean posterior probability.</p>
    </section>
    <section>
      <h2>Summary</h2>
      <div class="summary-grid">
        {"".join(f'<dl class="summary-card"><dt>{escape(label)}</dt><dd>{escape(value)}</dd></dl>' for label, value in summary_cards)}
      </div>
    </section>
    <section>
      <h2>Legend</h2>
      <div class="legend-bar"></div>
      <div class="legend-labels">
        <span>{escape(report.baseline_state or "0.0")}</span>
        <span>{escape(report.focal_state or "1.0")}</span>
      </div>
    </section>
    <section>
      <h2>Figure</h2>
      <div class="figure-shell">{svg_markup}</div>
      <p>Layout: <code>{escape(layout)}</code></p>
    </section>
  </main>
</body>
</html>
"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")


def render_stochastic_map_density_artifact(
    report: StochasticMapDensityReport,
    *,
    tree_path: Path,
    out_path: Path,
    layout: str = "phylogram",
) -> StochasticMapDensityArtifactResult:
    """Render one reviewer-facing stochastic-map density artifact."""
    if report.focal_state is None:
        raise ValueError(
            "stochastic-map density rendering requires one resolved focal state"
        )
    output_format = out_path.suffix.lower().lstrip(".")
    if output_format not in {"svg", "html"}:
        raise ValueError("stochastic-map density output must end in .svg or .html")
    branch_colors = {
        row.child_node: _blend_density_color(row.mean_posterior_probability)
        for row in report.branch_rows
    }
    svg_path = out_path if output_format == "svg" else out_path.with_suffix(".svg")
    render_result = render_tree_svg(
        tree_path,
        out_path=svg_path,
        layout=layout,
        branch_colors=branch_colors,
    )
    if output_format == "html":
        _write_stochastic_map_density_html(
            report=report,
            svg_path=svg_path,
            out_path=out_path,
            layout=layout,
        )
    return StochasticMapDensityArtifactResult(
        output_path=out_path,
        svg_path=svg_path,
        format=output_format,
        layout=layout,
        focal_state=report.focal_state,
        baseline_state=report.baseline_state,
        branch_count=len(report.branch_rows),
        rendered_branch_color_count=render_result.rendered_branch_color_count,
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
    mean_lookup = {row.transition: row.mean_count for row in report.aggregate_rows}
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


def write_stochastic_map_branch_probability_table(
    path: Path, report: StochasticMapDensityReport
) -> Path:
    """Export one per-branch state-probability table for a stochastic-map collection."""
    rows = [
        {
            "branch_index": row.branch_index,
            "parent_node": row.parent_node,
            "child_node": row.child_node,
            "state": row.state,
            "branch_length": row.branch_length,
            "mean_probability": row.mean_probability,
            "lower_95_probability": row.lower_95_probability,
            "upper_95_probability": row.upper_95_probability,
            "minimum_probability": row.minimum_probability,
            "maximum_probability": row.maximum_probability,
            "presence_fraction": row.presence_fraction,
        }
        for row in report.branch_state_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "branch_index",
            "parent_node",
            "child_node",
            "state",
            "branch_length",
            "mean_probability",
            "lower_95_probability",
            "upper_95_probability",
            "minimum_probability",
            "maximum_probability",
            "presence_fraction",
        ],
        rows=rows,
    )


def write_stochastic_map_density_branch_table(
    path: Path, report: StochasticMapDensityReport
) -> Path:
    """Export one per-branch focal-state density summary table."""
    rows = [
        {
            "branch_index": row.branch_index,
            "parent_node": row.parent_node,
            "child_node": row.child_node,
            "branch_length": row.branch_length,
            "focal_state": row.focal_state,
            "baseline_state": row.baseline_state or "",
            "mean_posterior_probability": row.mean_posterior_probability,
            "minimum_posterior_probability": row.minimum_posterior_probability,
            "maximum_posterior_probability": row.maximum_posterior_probability,
            "uncertainty": row.uncertainty,
            "slice_count": row.slice_count,
        }
        for row in report.branch_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "branch_index",
            "parent_node",
            "child_node",
            "branch_length",
            "focal_state",
            "baseline_state",
            "mean_posterior_probability",
            "minimum_posterior_probability",
            "maximum_posterior_probability",
            "uncertainty",
            "slice_count",
        ],
        rows=rows,
    )


def write_stochastic_map_density_slice_table(
    path: Path, report: StochasticMapDensityReport
) -> Path:
    """Export one flat branch-slice density table for a stochastic-map collection."""
    rows = [
        {
            "branch_index": row.branch_index,
            "parent_node": row.parent_node,
            "child_node": row.child_node,
            "branch_length": row.branch_length,
            "slice_index": row.slice_index,
            "start_depth": row.start_depth,
            "end_depth": row.end_depth,
            "start_time_fraction": row.start_time_fraction,
            "end_time_fraction": row.end_time_fraction,
            "posterior_probability": row.posterior_probability,
            "posterior_uncertainty": row.posterior_uncertainty,
        }
        for row in report.density_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "branch_index",
            "parent_node",
            "child_node",
            "branch_length",
            "slice_index",
            "start_depth",
            "end_depth",
            "start_time_fraction",
            "end_time_fraction",
            "posterior_probability",
            "posterior_uncertainty",
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
