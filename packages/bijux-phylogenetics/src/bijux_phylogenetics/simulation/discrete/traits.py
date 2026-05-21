from __future__ import annotations

from pathlib import Path
import random

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree


def _simulate_symmetric_state_trajectory(
    state: str,
    *,
    branch_length: float,
    rate: float,
    states: tuple[str, ...],
    rng: random.Random,
):
    from .._statistics import (
        _round_float,
    )
    from .._stochastic import _poisson_count
    from ..contracts import (
        SimulatedDiscreteStateSegment,
        SimulatedDiscreteTransitionEvent,
    )

    current_state = state
    events: list[SimulatedDiscreteTransitionEvent] = []
    event_count = _poisson_count(rate * branch_length, rng)
    event_times = (
        sorted(rng.random() * branch_length for _ in range(event_count))
        if event_count > 0
        else []
    )
    segments: list[SimulatedDiscreteStateSegment] = []
    previous_time = 0.0
    for index, event_time in enumerate(event_times, start=1):
        segments.append(
            SimulatedDiscreteStateSegment(
                parent_node="",
                child_node="",
                state=current_state,
                start_distance=_round_float(previous_time),
                end_distance=_round_float(event_time),
                duration=_round_float(event_time - previous_time),
            )
        )
        next_state = rng.choice(
            [candidate for candidate in states if candidate != current_state]
        )
        events.append(
            SimulatedDiscreteTransitionEvent(
                parent_node="",
                child_node="",
                source_state=current_state,
                target_state=next_state,
                event_index=index,
                branch_distance=_round_float(event_time),
            )
        )
        current_state = next_state
        previous_time = event_time
    segments.append(
        SimulatedDiscreteStateSegment(
            parent_node="",
            child_node="",
            state=current_state,
            start_distance=_round_float(previous_time),
            end_distance=_round_float(branch_length),
            duration=_round_float(branch_length - previous_time),
        )
    )
    return current_state, events, segments


def simulate_discrete_traits(
    tree_path: Path,
    *,
    states: list[str],
    transition_rate: float = 1.0,
    root_state: str | None = None,
    seed: int = 1,
):
    from bijux_phylogenetics.ancestral.common import (
        node_descendant_taxa,
        node_signature,
    )

    from .._state_propagation import _tip_values_from_node_map
    from ..contracts import (
        DiscreteHistoryRateRow,
        DiscreteTraitSimulationReport,
        SimulatedDiscreteBranchHistory,
        SimulatedDiscreteNode,
        SimulatedDiscreteStateSegment,
        SimulatedDiscreteTrait,
        SimulatedDiscreteTransitionEvent,
    )
    from .policy import (
        _normalize_discrete_states,
        _normalize_root_state_probabilities,
    )

    unique_states = _normalize_discrete_states(states)
    if transition_rate < 0.0:
        raise ValueError(f"transition_rate must be nonnegative, got {transition_rate}")
    tree = load_tree(tree_path)
    rng = random.Random(seed)  # nosec B311
    starting_state = root_state or unique_states[0]
    if starting_state not in unique_states:
        raise ValueError(f"root_state '{starting_state}' is not present in states")
    node_values: dict[str, str] = {}
    branch_histories: list[SimulatedDiscreteBranchHistory] = []

    def visit(node, state: str) -> None:
        parent_signature = node_signature(node)
        node_values[parent_signature] = state
        if node.is_leaf():
            return
        for child in node.children:
            branch_length = max(child.branch_length or 0.0, 0.0)
            child_signature = node_signature(child)
            child_state, events, segments = _simulate_symmetric_state_trajectory(
                state,
                branch_length=branch_length,
                rate=transition_rate,
                states=unique_states,
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
                    events=[
                        SimulatedDiscreteTransitionEvent(
                            parent_node=parent_signature,
                            child_node=child_signature,
                            source_state=event.source_state,
                            target_state=event.target_state,
                            event_index=event.event_index,
                            branch_distance=event.branch_distance,
                        )
                        for event in events
                    ],
                    segments=[
                        SimulatedDiscreteStateSegment(
                            parent_node=parent_signature,
                            child_node=child_signature,
                            state=segment.state,
                            start_distance=segment.start_distance,
                            end_distance=segment.end_distance,
                            duration=segment.duration,
                        )
                        for segment in segments
                    ],
                )
            )
            visit(child, child_state)

    visit(tree.root, starting_state)
    values = _tip_values_from_node_map(tree, node_values)
    return DiscreteTraitSimulationReport(
        model="symmetric-discrete",
        tree_path=tree_path,
        tip_count=tree.tip_count,
        seed=seed,
        states=list(unique_states),
        transition_rate=transition_rate,
        root_state=starting_state,
        root_state_probabilities=_normalize_root_state_probabilities(
            states=unique_states,
            root_state=starting_state,
            root_state_probabilities=None,
        ),
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
        rate_rows=[
            DiscreteHistoryRateRow(
                source_state=source_state,
                target_state=target_state,
                rate=float(transition_rate),
            )
            for source_state in unique_states
            for target_state in unique_states
            if source_state != target_state
        ],
    )


def write_discrete_trait_table(path: Path, report) -> Path:
    return write_taxon_rows(
        path,
        columns=["taxon", "state"],
        rows=[{"taxon": row.taxon, "state": row.state} for row in report.traits],
    )
