from __future__ import annotations

from dataclasses import dataclass, field
import math
from math import exp, sqrt
from pathlib import Path
import random
from statistics import median

from bijux_phylogenetics.core.alignment import AlignmentRecord
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.io.trees import load_tree


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

def _round_float(value: float) -> float:
    return round(float(value), 15)


def _mean(values: list[float]) -> float:
    return _round_float(sum(values) / len(values))


def _median(values: list[float]) -> float:
    return _round_float(float(median(values)))


def _population_standard_deviation(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mean_value = sum(values) / len(values)
    variance = sum((value - mean_value) ** 2 for value in values) / len(values)
    return _round_float(variance**0.5)


def _sample_standard_deviation(values: list[float]) -> float:
    if len(values) <= 1:
        return 0.0
    center = _mean(values)
    return _round_float(
        sqrt(sum((value - center) ** 2 for value in values) / (len(values) - 1))
    )


def _sample_covariance(left: list[float], right: list[float]) -> float:
    if len(left) != len(right):
        raise ValueError("covariance inputs must have the same length")
    if len(left) <= 1:
        return 0.0
    left_center = _mean(left)
    right_center = _mean(right)
    return _round_float(
        sum(
            (left_value - left_center) * (right_value - right_center)
            for left_value, right_value in zip(left, right, strict=True)
        )
        / (len(left) - 1)
    )


def _sample_correlation(left: list[float], right: list[float]) -> float:
    left_standard_deviation = _sample_standard_deviation(left)
    right_standard_deviation = _sample_standard_deviation(right)
    if left_standard_deviation == 0.0 or right_standard_deviation == 0.0:
        return 0.0
    return _round_float(
        _sample_covariance(left, right)
        / (left_standard_deviation * right_standard_deviation)
    )


def simulate_birth_death_trees(
    *,
    tree_count: int,
    tip_count: int,
    birth_rate: float = 1.0,
    death_rate: float = 0.25,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
) -> tuple[list[PhyloTree], TreeSimulationReport]:
    from .trees import simulate_birth_death_trees as simulate_birth_death_trees_impl

    return simulate_birth_death_trees_impl(
        tree_count=tree_count,
        tip_count=tip_count,
        birth_rate=birth_rate,
        death_rate=death_rate,
        seed=seed,
        taxon_prefix=taxon_prefix,
    )


def simulate_random_trees(
    *,
    tree_count: int,
    tip_count: int,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
    branch_length_model: str = "uniform",
) -> tuple[list[PhyloTree], TreeSimulationReport]:
    from .trees import simulate_random_trees as simulate_random_trees_impl

    return simulate_random_trees_impl(
        tree_count=tree_count,
        tip_count=tip_count,
        seed=seed,
        taxon_prefix=taxon_prefix,
        branch_length_model=branch_length_model,
    )


def simulate_random_tree(
    *,
    tip_count: int,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
    branch_length_model: str = "uniform",
) -> tuple[PhyloTree, TreeSimulationReport]:
    from .trees import simulate_random_tree as simulate_random_tree_impl

    return simulate_random_tree_impl(
        tip_count=tip_count,
        seed=seed,
        taxon_prefix=taxon_prefix,
        branch_length_model=branch_length_model,
    )




def _iter_tip_trait_values(
    tree: PhyloTree,
    *,
    root_state: float,
    propagate,
) -> dict[str, float]:
    values: dict[str, float] = {}

    def visit(node: TreeNode, state: float) -> None:
        if node.is_leaf():
            if node.name is not None:
                values[node.name] = (
                    round(state, 15) if isinstance(state, float) else state
                )
            return
        for child in node.children:
            branch_length = max(child.branch_length or 0.0, 0.0)
            visit(child, propagate(state, branch_length))

    visit(tree.root, root_state)
    return values


def _iter_node_trait_values(
    tree: PhyloTree,
    *,
    root_state,
    propagate,
) -> dict[str, object]:
    from bijux_phylogenetics.ancestral.common import node_signature

    values: dict[str, object] = {}

    def visit(node: TreeNode, state) -> None:
        values[node_signature(node)] = state
        if node.is_leaf():
            return
        for child in node.children:
            branch_length = max(child.branch_length or 0.0, 0.0)
            visit(child, propagate(state, branch_length))

    visit(tree.root, root_state)
    return values


def _tip_values_from_node_map(
    tree: PhyloTree, node_values: dict[str, object]
) -> dict[str, object]:
    from bijux_phylogenetics.ancestral.common import node_signature

    return {
        node.name: (
            round(float(node_values[node_signature(node)]), 15)
            if isinstance(node_values[node_signature(node)], float)
            else node_values[node_signature(node)]
        )
        for node in tree.iter_leaves()
        if node.name is not None
    }


def _resolve_brownian_sigma_parameters(
    *,
    sigma: float | None,
    sigma_squared: float | None,
) -> tuple[float, float]:
    if sigma is None and sigma_squared is None:
        return 1.0, 1.0
    if sigma is not None and sigma < 0.0:
        raise ValueError(f"sigma must be nonnegative, got {sigma}")
    if sigma_squared is not None and sigma_squared < 0.0:
        raise ValueError(f"sigma_squared must be nonnegative, got {sigma_squared}")
    if sigma is None:
        resolved_sigma_squared = float(sigma_squared)
        return sqrt(resolved_sigma_squared), resolved_sigma_squared
    resolved_sigma_squared = sigma * sigma
    if sigma_squared is None:
        return sigma, resolved_sigma_squared
    if not math.isclose(
        resolved_sigma_squared,
        sigma_squared,
        rel_tol=1e-12,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "sigma and sigma_squared must describe the same Brownian rate parameter"
        )
    return sigma, float(sigma_squared)


def _simulate_brownian_node_values(
    tree: PhyloTree,
    *,
    root_state: float,
    sigma: float,
    rng: random.Random,
) -> dict[str, float]:
    return _iter_node_trait_values(
        tree,
        root_state=root_state,
        propagate=lambda state, branch_length: (
            state + rng.gauss(0.0, sigma * sqrt(branch_length))
        ),
    )


def _poisson_count(expected_changes: float, rng: random.Random) -> int:
    if expected_changes <= 0.0:
        return 0
    threshold = exp(-expected_changes)
    product = 1.0
    changes = 0
    while product > threshold:
        changes += 1
        product *= rng.random()
    return changes - 1


def _simulate_symmetric_state_trajectory(
    state: str,
    *,
    branch_length: float,
    rate: float,
    states: tuple[str, ...],
    rng: random.Random,
) -> tuple[
    str,
    list[SimulatedDiscreteTransitionEvent],
    list[SimulatedDiscreteStateSegment],
]:
    if rate < 0.0:
        raise ValueError(f"rate must be nonnegative, got {rate}")
    next_state = state
    event_states: list[tuple[str, str]] = []
    segment_boundaries = [0.0]
    for _ in range(_poisson_count(rate * branch_length, rng)):
        alternatives = [candidate for candidate in states if candidate != next_state]
        candidate = rng.choice(alternatives)
        event_states.append((next_state, candidate))
        next_state = candidate
    if event_states:
        event_distances = [
            branch_length * index / (len(event_states) + 1)
            for index in range(1, len(event_states) + 1)
        ]
        segment_boundaries.extend(event_distances)
    segment_boundaries.append(branch_length)
    events: list[SimulatedDiscreteTransitionEvent] = []
    segments: list[SimulatedDiscreteStateSegment] = []
    current_state = state
    parent_node = ""
    child_node = ""
    for index, (start_distance, end_distance) in enumerate(
        zip(segment_boundaries[:-1], segment_boundaries[1:], strict=True),
        start=1,
    ):
        segments.append(
            SimulatedDiscreteStateSegment(
                parent_node=parent_node,
                child_node=child_node,
                state=current_state,
                start_distance=_round_float(start_distance),
                end_distance=_round_float(end_distance),
                duration=_round_float(max(end_distance - start_distance, 0.0)),
            )
        )
        if index <= len(event_states):
            source_state, target_state = event_states[index - 1]
            events.append(
                SimulatedDiscreteTransitionEvent(
                    parent_node=parent_node,
                    child_node=child_node,
                    source_state=source_state,
                    target_state=target_state,
                    event_index=index,
                    branch_distance=_round_float(end_distance),
                )
            )
            current_state = target_state
    return next_state, events, segments


def _quantile(values: list[float], probability: float) -> float:
    if not values:
        return 0.0
    if len(values) == 1:
        return _round_float(values[0])
    sorted_values = sorted(values)
    scaled_index = (len(sorted_values) - 1) * probability
    lower_index = int(scaled_index)
    upper_index = min(lower_index + 1, len(sorted_values) - 1)
    fraction = scaled_index - lower_index
    lower_value = sorted_values[lower_index]
    upper_value = sorted_values[upper_index]
    return _round_float(lower_value + (upper_value - lower_value) * fraction)


def _normalize_discrete_states(states: list[str]) -> tuple[str, ...]:
    unique_states = tuple(dict.fromkeys(state for state in states if state))
    if len(unique_states) < 2:
        raise ValueError("states must contain at least two distinct non-empty states")
    return unique_states


def _normalize_rate_rows(
    *,
    states: tuple[str, ...],
    rate_rows: list[DiscreteHistoryRateRow],
) -> list[DiscreteHistoryRateRow]:
    state_set = set(states)
    normalized_rows: list[DiscreteHistoryRateRow] = []
    seen_pairs: set[tuple[str, str]] = set()
    for row in rate_rows:
        if row.source_state not in state_set:
            raise ValueError(
                f"unknown source_state '{row.source_state}' in rate matrix"
            )
        if row.target_state not in state_set:
            raise ValueError(
                f"unknown target_state '{row.target_state}' in rate matrix"
            )
        if row.source_state == row.target_state:
            raise ValueError("rate matrix rows must describe source->target changes")
        if row.rate < 0.0:
            raise ValueError(f"rate must be nonnegative, got {row.rate}")
        pair = (row.source_state, row.target_state)
        if pair in seen_pairs:
            raise ValueError(
                f"duplicate rate row for {row.source_state}->{row.target_state}"
            )
        seen_pairs.add(pair)
        normalized_rows.append(
            DiscreteHistoryRateRow(
                source_state=row.source_state,
                target_state=row.target_state,
                rate=_round_float(row.rate),
            )
        )
    return sorted(normalized_rows, key=lambda row: (row.source_state, row.target_state))


def _normalize_root_state_probabilities(
    *,
    states: tuple[str, ...],
    root_state: str | None,
    root_state_probabilities: dict[str, float] | None,
) -> dict[str, float]:
    if root_state is not None and root_state_probabilities is not None:
        raise ValueError(
            "root_state and root_state_probabilities cannot be supplied together"
        )
    state_set = set(states)
    if root_state is not None:
        if root_state not in state_set:
            raise ValueError(f"root_state '{root_state}' is not present in states")
        return {state: 1.0 if state == root_state else 0.0 for state in states}
    if root_state_probabilities is None:
        probability = 1.0 / len(states)
        return {state: _round_float(probability) for state in states}
    unknown_states = set(root_state_probabilities).difference(state_set)
    if unknown_states:
        unknown_state = sorted(unknown_states)[0]
        raise ValueError(
            f"root_state_probabilities contains unknown state '{unknown_state}'"
        )
    probabilities = {
        state: float(root_state_probabilities.get(state, 0.0)) for state in states
    }
    if any(value < 0.0 for value in probabilities.values()):
        raise ValueError("root_state_probabilities cannot contain negative values")
    total = sum(probabilities.values())
    if total <= 0.0:
        raise ValueError("root_state_probabilities must sum to a positive value")
    return {
        state: _round_float(value / total) for state, value in probabilities.items()
    }


def _sample_discrete_state(
    probabilities: dict[str, float],
    *,
    rng: random.Random,
) -> str:
    threshold = rng.random()
    cumulative = 0.0
    last_state = next(iter(probabilities))
    for state, probability in probabilities.items():
        cumulative += probability
        last_state = state
        if threshold <= cumulative:
            return state
    return last_state


def _build_rate_lookup(
    *,
    states: tuple[str, ...],
    rate_rows: list[DiscreteHistoryRateRow],
) -> dict[str, dict[str, float]]:
    lookup = {
        source_state: {
            target_state: 0.0 for target_state in states if target_state != source_state
        }
        for source_state in states
    }
    for row in rate_rows:
        lookup[row.source_state][row.target_state] = row.rate
    return lookup


def _simulate_rate_matrix_state_trajectory(
    state: str,
    *,
    parent_node: str,
    child_node: str,
    branch_length: float,
    rate_lookup: dict[str, dict[str, float]],
    rng: random.Random,
) -> tuple[
    str,
    list[SimulatedDiscreteTransitionEvent],
    list[SimulatedDiscreteStateSegment],
]:
    current_state = state
    elapsed = 0.0
    events: list[SimulatedDiscreteTransitionEvent] = []
    segments: list[SimulatedDiscreteStateSegment] = []
    while elapsed < branch_length:
        outgoing = rate_lookup[current_state]
        exit_rate = sum(outgoing.values())
        if exit_rate <= 0.0:
            segments.append(
                SimulatedDiscreteStateSegment(
                    parent_node=parent_node,
                    child_node=child_node,
                    state=current_state,
                    start_distance=_round_float(elapsed),
                    end_distance=_round_float(branch_length),
                    duration=_round_float(branch_length - elapsed),
                )
            )
            elapsed = branch_length
            break
        wait_time = rng.expovariate(exit_rate)
        next_time = min(branch_length, elapsed + wait_time)
        segments.append(
            SimulatedDiscreteStateSegment(
                parent_node=parent_node,
                child_node=child_node,
                state=current_state,
                start_distance=_round_float(elapsed),
                end_distance=_round_float(next_time),
                duration=_round_float(max(next_time - elapsed, 0.0)),
            )
        )
        if next_time >= branch_length:
            elapsed = branch_length
            break
        next_state = _sample_discrete_state(
            {
                target_state: rate / exit_rate
                for target_state, rate in outgoing.items()
                if rate > 0.0
            },
            rng=rng,
        )
        events.append(
            SimulatedDiscreteTransitionEvent(
                parent_node=parent_node,
                child_node=child_node,
                source_state=current_state,
                target_state=next_state,
                event_index=len(events) + 1,
                branch_distance=_round_float(next_time),
            )
        )
        current_state = next_state
        elapsed = next_time
    if not segments:
        segments.append(
            SimulatedDiscreteStateSegment(
                parent_node=parent_node,
                child_node=child_node,
                state=current_state,
                start_distance=0.0,
                end_distance=_round_float(branch_length),
                duration=_round_float(branch_length),
            )
        )
    return current_state, events, segments


def _simulate_discrete_history_once(
    tree_path: Path,
    *,
    model: str,
    states: tuple[str, ...],
    rate_rows: list[DiscreteHistoryRateRow],
    transition_rate: float | None,
    fixed_root_state: str | None,
    root_state_probabilities: dict[str, float],
    seed: int,
) -> DiscreteTraitSimulationReport:
    from bijux_phylogenetics.ancestral.common import (
        node_descendant_taxa,
        node_signature,
    )

    tree = load_tree(tree_path)
    rng = random.Random(seed)  # nosec B311
    node_values: dict[str, str] = {}
    branch_histories: list[SimulatedDiscreteBranchHistory] = []
    rate_lookup = _build_rate_lookup(states=states, rate_rows=rate_rows)
    starting_state = (
        fixed_root_state
        if fixed_root_state is not None
        else _sample_discrete_state(root_state_probabilities, rng=rng)
    )

    def visit(node: TreeNode, state: str) -> None:
        parent_signature = node_signature(node)
        node_values[parent_signature] = state
        if node.is_leaf():
            return
        for child in node.children:
            branch_length = max(child.branch_length or 0.0, 0.0)
            child_signature = node_signature(child)
            child_state, events, segments = _simulate_rate_matrix_state_trajectory(
                state,
                parent_node=parent_signature,
                child_node=child_signature,
                branch_length=branch_length,
                rate_lookup=rate_lookup,
                rng=rng,
            )
            branch_histories.append(
                SimulatedDiscreteBranchHistory(
                    parent_node=parent_signature,
                    child_node=child_signature,
                    branch_length=float(format(branch_length, ".15g")),
                    start_state=state,
                    end_state=child_state,
                    changed=bool(events),
                    event_count=len(events),
                    events=events,
                    segments=segments,
                )
            )
            visit(child, child_state)

    visit(tree.root, starting_state)
    values = _tip_values_from_node_map(tree, node_values)
    return DiscreteTraitSimulationReport(
        model=model,
        tree_path=tree_path,
        tip_count=tree.tip_count,
        seed=seed,
        states=list(states),
        transition_rate=transition_rate,
        root_state=starting_state,
        root_state_probabilities=dict(root_state_probabilities),
        traits=[
            SimulatedDiscreteTrait(taxon=taxon, state=state)
            for taxon, state in sorted(values.items())
        ],
        node_states=[
            SimulatedDiscreteNode(
                node=node_signature(node),
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=node_descendant_taxa(node),
                state=str(node_values[node_signature(node)]),
            )
            for node in tree.iter_nodes()
        ],
        branch_histories=branch_histories,
        rate_rows=list(rate_rows),
    )


def _summarize_discrete_history_collection(
    simulations: list[DiscreteTraitSimulationReport],
    *,
    states: tuple[str, ...],
) -> tuple[float, float, float, list[DiscreteHistorySummaryRow]]:
    total_transition_counts = [
        float(sum(branch.event_count for branch in simulation.branch_histories))
        for simulation in simulations
    ]
    transition_labels = [
        f"{source_state}->{target_state}"
        for source_state in states
        for target_state in states
        if source_state != target_state
    ]
    transition_values = {label: [] for label in transition_labels}
    state_time_values = {state: [] for state in states}
    tip_state_values: dict[str, list[float]] = {}
    for simulation in simulations:
        branch_transition_counts = dict.fromkeys(transition_labels, 0.0)
        branch_state_totals = dict.fromkeys(states, 0.0)
        tip_lookup = {row.taxon: row.state for row in simulation.traits}
        for history in simulation.branch_histories:
            for event in history.events:
                branch_transition_counts[
                    f"{event.source_state}->{event.target_state}"
                ] += 1.0
            for segment in history.segments:
                branch_state_totals[segment.state] += segment.duration
        for label, value in branch_transition_counts.items():
            transition_values[label].append(value)
        for state, value in branch_state_totals.items():
            state_time_values[state].append(value)
        for taxon, state in sorted(tip_lookup.items()):
            for candidate_state in states:
                key = f"{taxon}:{candidate_state}"
                tip_state_values.setdefault(key, []).append(
                    1.0 if state == candidate_state else 0.0
                )
    rows: list[DiscreteHistorySummaryRow] = []
    for label, values in sorted(transition_values.items()):
        rows.append(
            DiscreteHistorySummaryRow(
                row_kind="transition_count",
                label=label,
                mean_value=_mean(values),
                lower_95_interval=_quantile(values, 0.025),
                upper_95_interval=_quantile(values, 0.975),
                presence_fraction=_round_float(
                    sum(1 for value in values if value > 0.0) / max(len(values), 1)
                ),
            )
        )
    for state, values in sorted(state_time_values.items()):
        rows.append(
            DiscreteHistorySummaryRow(
                row_kind="state_time",
                label=state,
                mean_value=_mean(values),
                lower_95_interval=_quantile(values, 0.025),
                upper_95_interval=_quantile(values, 0.975),
                # Every replicate has a defined total time for every state, even when that
                # total is zero, so the row is always present in the summary surface.
                presence_fraction=1.0,
            )
        )
    for label, values in sorted(tip_state_values.items()):
        rows.append(
            DiscreteHistorySummaryRow(
                row_kind="tip_state_frequency",
                label=label,
                mean_value=_mean(values),
                lower_95_interval=_quantile(values, 0.025),
                upper_95_interval=_quantile(values, 0.975),
                presence_fraction=_mean(values),
            )
        )
    return (
        _mean(total_transition_counts),
        _quantile(total_transition_counts, 0.025),
        _quantile(total_transition_counts, 0.975),
        rows,
    )


def simulate_discrete_histories(
    tree_path: Path,
    *,
    states: list[str],
    rate_rows: list[DiscreteHistoryRateRow],
    root_state: str | None = None,
    root_state_probabilities: dict[str, float] | None = None,
    replicates: int = 1,
    seed: int = 1,
) -> DiscreteHistorySimulationCollectionReport:
    """Simulate one or more discrete histories on a fixed tree from one rate matrix."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    normalized_states = _normalize_discrete_states(states)
    normalized_rates = _normalize_rate_rows(
        states=normalized_states,
        rate_rows=rate_rows,
    )
    normalized_root_probabilities = _normalize_root_state_probabilities(
        states=normalized_states,
        root_state=root_state,
        root_state_probabilities=root_state_probabilities,
    )
    simulations = [
        _simulate_discrete_history_once(
            tree_path,
            model="rate-matrix-discrete-history",
            states=normalized_states,
            rate_rows=normalized_rates,
            transition_rate=None,
            fixed_root_state=root_state,
            root_state_probabilities=normalized_root_probabilities,
            seed=seed + index - 1,
        )
        for index in range(1, replicates + 1)
    ]
    tree = load_tree(tree_path)
    (
        mean_total_transition_count,
        lower_95_total_transition_count,
        upper_95_total_transition_count,
        rows,
    ) = _summarize_discrete_history_collection(
        simulations,
        states=normalized_states,
    )
    return DiscreteHistorySimulationCollectionReport(
        model="rate-matrix-discrete-history",
        tree_path=tree_path,
        tip_count=tree.tip_count,
        branch_count=sum(1 for _ in tree.iter_nodes()) - 1,
        replicate_count=replicates,
        seed=seed,
        states=list(normalized_states),
        fixed_root_state=root_state,
        root_state_probabilities=normalized_root_probabilities,
        rate_rows=normalized_rates,
        simulations=simulations,
        mean_total_transition_count=mean_total_transition_count,
        lower_95_total_transition_count=lower_95_total_transition_count,
        upper_95_total_transition_count=upper_95_total_transition_count,
        rows=rows,
    )


def simulate_coalescent_trees(
    *,
    tree_count: int,
    tip_count: int,
    population_size: float = 1.0,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
) -> tuple[list[PhyloTree], TreeSimulationReport]:
    from .trees import simulate_coalescent_trees as simulate_coalescent_trees_impl

    return simulate_coalescent_trees_impl(
        tree_count=tree_count,
        tip_count=tip_count,
        population_size=population_size,
        seed=seed,
        taxon_prefix=taxon_prefix,
    )


def simulate_coalescent_tree(
    *,
    tip_count: int,
    population_size: float = 1.0,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
) -> tuple[PhyloTree, TreeSimulationReport]:
    from .trees import simulate_coalescent_tree as simulate_coalescent_tree_impl

    return simulate_coalescent_tree_impl(
        tip_count=tip_count,
        population_size=population_size,
        seed=seed,
        taxon_prefix=taxon_prefix,
    )


def simulate_brownian_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float | None = None,
    sigma_squared: float | None = None,
    seed: int = 1,
) -> ContinuousTraitSimulationReport:
    from .continuous import simulate_brownian_traits as simulate_brownian_traits_impl

    return simulate_brownian_traits_impl(
        tree_path,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma_squared,
        seed=seed,
    )


def simulate_ou_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float = 1.0,
    alpha: float = 1.0,
    theta: float = 0.0,
    seed: int = 1,
) -> ContinuousTraitSimulationReport:
    from .continuous import simulate_ou_traits as simulate_ou_traits_impl

    return simulate_ou_traits_impl(
        tree_path,
        root_state=root_state,
        sigma=sigma,
        alpha=alpha,
        theta=theta,
        seed=seed,
    )


def simulate_early_burst_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float = 1.0,
    rate_change: float = 1.0,
    seed: int = 1,
) -> ContinuousTraitSimulationReport:
    from .continuous import (
        simulate_early_burst_traits as simulate_early_burst_traits_impl,
    )

    return simulate_early_burst_traits_impl(
        tree_path,
        root_state=root_state,
        sigma=sigma,
        rate_change=rate_change,
        seed=seed,
    )


def simulate_brownian_trait_collection(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float | None = None,
    sigma_squared: float | None = None,
    replicates: int = 128,
    seed: int = 1,
) -> ContinuousTraitSimulationCollectionReport:
    from .continuous import (
        simulate_brownian_trait_collection as simulate_brownian_trait_collection_impl,
    )

    return simulate_brownian_trait_collection_impl(
        tree_path,
        root_state=root_state,
        sigma=sigma,
        sigma_squared=sigma_squared,
        replicates=replicates,
        seed=seed,
    )


def simulate_correlated_brownian_traits(
    tree_path: Path,
    *,
    trait_names: list[str] | tuple[str, ...],
    evolutionary_covariance_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    evolutionary_correlation_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    trait_standard_deviations: list[float] | tuple[float, ...] | None = None,
    root_states: list[float] | tuple[float, ...] | None = None,
    seed: int = 1,
) -> CorrelatedContinuousTraitSimulationReport:
    from .correlated import (
        simulate_correlated_brownian_traits as simulate_correlated_brownian_traits_impl,
    )

    return simulate_correlated_brownian_traits_impl(
        tree_path,
        trait_names=trait_names,
        evolutionary_covariance_matrix=evolutionary_covariance_matrix,
        evolutionary_correlation_matrix=evolutionary_correlation_matrix,
        trait_standard_deviations=trait_standard_deviations,
        root_states=root_states,
        seed=seed,
    )


def simulate_correlated_brownian_trait_collection(
    tree_path: Path,
    *,
    trait_names: list[str] | tuple[str, ...],
    evolutionary_covariance_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    evolutionary_correlation_matrix: list[list[float]]
    | tuple[tuple[float, ...], ...]
    | None = None,
    trait_standard_deviations: list[float] | tuple[float, ...] | None = None,
    root_states: list[float] | tuple[float, ...] | None = None,
    replicates: int = 128,
    seed: int = 1,
) -> CorrelatedContinuousTraitSimulationCollectionReport:
    from .correlated import (
        simulate_correlated_brownian_trait_collection as simulate_correlated_brownian_trait_collection_impl,
    )

    return simulate_correlated_brownian_trait_collection_impl(
        tree_path,
        trait_names=trait_names,
        evolutionary_covariance_matrix=evolutionary_covariance_matrix,
        evolutionary_correlation_matrix=evolutionary_correlation_matrix,
        trait_standard_deviations=trait_standard_deviations,
        root_states=root_states,
        replicates=replicates,
        seed=seed,
    )


def simulate_discrete_traits(
    tree_path: Path,
    *,
    states: list[str],
    transition_rate: float = 1.0,
    root_state: str | None = None,
    seed: int = 1,
) -> DiscreteTraitSimulationReport:
    from .discrete_traits import (
        simulate_discrete_traits as simulate_discrete_traits_impl,
    )

    return simulate_discrete_traits_impl(
        tree_path,
        states=states,
        transition_rate=transition_rate,
        root_state=root_state,
        seed=seed,
    )


def simulate_dna_alignment(
    tree_path: Path,
    *,
    sequence_length: int,
    substitution_rate: float = 1.0,
    seed: int = 1,
) -> AlignmentSimulationReport:
    from .alignment import simulate_dna_alignment as simulate_dna_alignment_impl

    return simulate_dna_alignment_impl(
        tree_path,
        sequence_length=sequence_length,
        substitution_rate=substitution_rate,
        seed=seed,
    )


def simulate_protein_alignment(
    tree_path: Path,
    *,
    sequence_length: int,
    substitution_rate: float = 1.0,
    seed: int = 1,
) -> AlignmentSimulationReport:
    from .alignment import (
        simulate_protein_alignment as simulate_protein_alignment_impl,
    )

    return simulate_protein_alignment_impl(
        tree_path,
        sequence_length=sequence_length,
        substitution_rate=substitution_rate,
        seed=seed,
    )


def write_tree_set(path: Path, trees: list[PhyloTree]) -> Path:
    from .trees import write_tree_set as write_tree_set_impl

    return write_tree_set_impl(path, trees)


def write_simulated_tree(path: Path, tree: PhyloTree) -> Path:
    from .trees import write_simulated_tree as write_simulated_tree_impl

    return write_simulated_tree_impl(path, tree)


def write_tree_simulation_record_table(
    path: Path, report: TreeSimulationReport
) -> Path:
    from .trees import (
        write_tree_simulation_record_table as write_tree_simulation_record_table_impl,
    )

    return write_tree_simulation_record_table_impl(path, report)


def write_tree_simulation_envelope_table(
    path: Path,
    report: TreeSimulationReport,
) -> Path:
    from .trees import (
        write_tree_simulation_envelope_table as write_tree_simulation_envelope_table_impl,
    )

    return write_tree_simulation_envelope_table_impl(path, report)


def write_continuous_trait_table(
    path: Path, report: ContinuousTraitSimulationReport
) -> Path:
    from .continuous import write_continuous_trait_table as write_continuous_trait_table_impl

    return write_continuous_trait_table_impl(path, report)


def write_continuous_trait_collection_table(
    path: Path,
    report: ContinuousTraitSimulationCollectionReport,
) -> Path:
    from .continuous import (
        write_continuous_trait_collection_table as write_continuous_trait_collection_table_impl,
    )

    return write_continuous_trait_collection_table_impl(path, report)


def write_continuous_trait_collection_summary_table(
    path: Path,
    report: ContinuousTraitSimulationCollectionReport,
) -> Path:
    from .continuous import (
        write_continuous_trait_collection_summary_table as write_continuous_trait_collection_summary_table_impl,
    )

    return write_continuous_trait_collection_summary_table_impl(path, report)


def write_correlated_continuous_trait_table(
    path: Path,
    report: CorrelatedContinuousTraitSimulationReport,
) -> Path:
    from .correlated import (
        write_correlated_continuous_trait_table as write_correlated_continuous_trait_table_impl,
    )

    return write_correlated_continuous_trait_table_impl(path, report)


def write_correlated_continuous_trait_collection_table(
    path: Path,
    report: CorrelatedContinuousTraitSimulationCollectionReport,
) -> Path:
    from .correlated import (
        write_correlated_continuous_trait_collection_table as write_correlated_continuous_trait_collection_table_impl,
    )

    return write_correlated_continuous_trait_collection_table_impl(path, report)


def write_correlated_continuous_trait_collection_summary_table(
    path: Path,
    report: CorrelatedContinuousTraitSimulationCollectionReport,
) -> Path:
    from .correlated import (
        write_correlated_continuous_trait_collection_summary_table as write_correlated_continuous_trait_collection_summary_table_impl,
    )

    return write_correlated_continuous_trait_collection_summary_table_impl(path, report)


def write_discrete_trait_table(
    path: Path, report: DiscreteTraitSimulationReport
) -> Path:
    from .discrete_traits import (
        write_discrete_trait_table as write_discrete_trait_table_impl,
    )

    return write_discrete_trait_table_impl(path, report)


def write_discrete_history_tip_truth_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    """Write simulated discrete tip states across history replicates."""
    return write_taxon_rows(
        path,
        columns=["replicate_index", "taxon", "state"],
        rows=[
            {
                "replicate_index": replicate_index,
                "taxon": row.taxon,
                "state": row.state,
            }
            for replicate_index, simulation in enumerate(report.simulations, start=1)
            for row in simulation.traits
        ],
    )


def write_discrete_history_node_truth_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    """Write simulated node states across history replicates."""
    return write_taxon_rows(
        path,
        columns=[
            "replicate_index",
            "node",
            "node_name",
            "is_tip",
            "descendant_taxa",
            "state",
        ],
        rows=[
            {
                "replicate_index": replicate_index,
                "node": row.node,
                "node_name": row.node_name or "",
                "is_tip": str(row.is_tip).lower(),
                "descendant_taxa": ",".join(row.descendant_taxa),
                "state": row.state,
            }
            for replicate_index, simulation in enumerate(report.simulations, start=1)
            for row in simulation.node_states
        ],
    )


def write_discrete_history_branch_truth_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    """Write one branch-history truth table across simulated discrete histories."""
    return write_taxon_rows(
        path,
        columns=[
            "replicate_index",
            "parent_node",
            "child_node",
            "branch_length",
            "start_state",
            "end_state",
            "changed",
            "event_count",
        ],
        rows=[
            {
                "replicate_index": replicate_index,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "branch_length": row.branch_length,
                "start_state": row.start_state,
                "end_state": row.end_state,
                "changed": str(row.changed).lower(),
                "event_count": row.event_count,
            }
            for replicate_index, simulation in enumerate(report.simulations, start=1)
            for row in simulation.branch_histories
        ],
    )


def write_discrete_history_event_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    """Write one flat discrete-history transition event ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "replicate_index",
            "parent_node",
            "child_node",
            "source_state",
            "target_state",
            "event_index",
            "branch_distance",
        ],
        rows=[
            {
                "replicate_index": replicate_index,
                "parent_node": branch.parent_node,
                "child_node": branch.child_node,
                "source_state": event.source_state,
                "target_state": event.target_state,
                "event_index": event.event_index,
                "branch_distance": event.branch_distance,
            }
            for replicate_index, simulation in enumerate(report.simulations, start=1)
            for branch in simulation.branch_histories
            for event in branch.events
        ],
    )


def write_discrete_history_segment_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    """Write one flat discrete-history branch-segment truth ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "replicate_index",
            "parent_node",
            "child_node",
            "state",
            "start_distance",
            "end_distance",
            "duration",
        ],
        rows=[
            {
                "replicate_index": replicate_index,
                "parent_node": branch.parent_node,
                "child_node": branch.child_node,
                "state": segment.state,
                "start_distance": segment.start_distance,
                "end_distance": segment.end_distance,
                "duration": segment.duration,
            }
            for replicate_index, simulation in enumerate(report.simulations, start=1)
            for branch in simulation.branch_histories
            for segment in branch.segments
        ],
    )


def write_discrete_history_summary_table(
    path: Path,
    report: DiscreteHistorySimulationCollectionReport,
) -> Path:
    """Write one discrete-history summary table for parity and recovery review."""
    return write_taxon_rows(
        path,
        columns=[
            "row_kind",
            "label",
            "mean_value",
            "lower_95_interval",
            "upper_95_interval",
            "presence_fraction",
        ],
        rows=[
            {
                "row_kind": row.row_kind,
                "label": row.label,
                "mean_value": row.mean_value,
                "lower_95_interval": row.lower_95_interval,
                "upper_95_interval": row.upper_95_interval,
                "presence_fraction": row.presence_fraction,
            }
            for row in report.rows
        ],
    )


def write_simulated_alignment(path: Path, report: AlignmentSimulationReport) -> Path:
    from .alignment import write_simulated_alignment as write_simulated_alignment_impl

    return write_simulated_alignment_impl(path, report)
