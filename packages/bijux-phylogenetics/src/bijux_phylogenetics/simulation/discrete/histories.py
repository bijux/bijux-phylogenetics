from __future__ import annotations

from pathlib import Path
import random
from typing import TYPE_CHECKING

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

if TYPE_CHECKING:
    from ..contracts import DiscreteHistoryRateRow, DiscreteTraitSimulationReport


def _quantile(values: list[float], probability: float) -> float:
    from .._statistics import _round_float

    if not values:
        return 0.0
    if len(values) == 1:
        return _round_float(values[0])
    ordered_values = sorted(values)
    position = (len(ordered_values) - 1) * probability
    lower_index = int(position)
    upper_index = min(lower_index + 1, len(ordered_values) - 1)
    fraction = position - lower_index
    interpolated = (
        ordered_values[lower_index]
        + (ordered_values[upper_index] - ordered_values[lower_index]) * fraction
    )
    return _round_float(interpolated)


def _normalize_rate_rows(
    *, states: tuple[str, ...], rate_rows: list[DiscreteHistoryRateRow]
) -> list[DiscreteHistoryRateRow]:
    from .._statistics import _round_float
    from ..contracts import DiscreteHistoryRateRow

    if not rate_rows:
        raise ValueError("rate_rows must contain at least one transition rate row")
    normalized_rows: list[DiscreteHistoryRateRow] = []
    seen_pairs: set[tuple[str, str]] = set()
    state_set = set(states)
    for row in rate_rows:
        if row.source_state not in state_set or row.target_state not in state_set:
            raise ValueError("rate matrix rows must use only declared states")
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
):
    from .._statistics import (
        _round_float,
    )
    from ..contracts import (
        SimulatedDiscreteStateSegment,
        SimulatedDiscreteTransitionEvent,
    )

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
                end_distance=branch_length,
                duration=branch_length,
            )
        )
    return current_state, events, segments


def _simulate_discrete_history_once(
    tree: PhyloTree,
    *,
    tree_path: Path,
    model: str,
    states: tuple[str, ...],
    rate_rows: list[DiscreteHistoryRateRow],
    transition_rate: float | None,
    fixed_root_state: str | None,
    root_state_probabilities: dict[str, float],
    transform_name: str | None,
    transform_parameter_name: str | None,
    transform_parameter_value: float | None,
    seed: int,
):
    from bijux_phylogenetics.ancestral.common import (
        node_descendant_taxa,
        node_signature,
    )

    from .._state_propagation import _tip_values_from_node_map
    from ..contracts import (
        DiscreteTraitSimulationReport,
        SimulatedDiscreteBranchHistory,
        SimulatedDiscreteNode,
        SimulatedDiscreteTrait,
    )

    rng = random.Random(seed)  # nosec B311
    node_values: dict[str, str] = {}
    branch_histories: list[SimulatedDiscreteBranchHistory] = []
    rate_lookup = _build_rate_lookup(states=states, rate_rows=rate_rows)
    starting_state = (
        fixed_root_state
        if fixed_root_state is not None
        else _sample_discrete_state(root_state_probabilities, rng=rng)
    )

    def visit(node, state: str) -> None:
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
        transform_name=transform_name,
        transform_parameter_name=transform_parameter_name,
        transform_parameter_value=transform_parameter_value,
    )


def _resolve_discrete_history_transform_name(transform: str | None) -> str | None:
    if transform is None:
        return None
    aliases = {
        "lambda": "lambda",
        "pagel-lambda": "lambda",
        "kappa": "kappa",
        "pagel-kappa": "kappa",
        "delta": "delta",
        "pagel-delta": "delta",
        "EB": "early-burst",
        "early-burst": "early-burst",
    }
    resolved = aliases.get(transform)
    if resolved is None:
        raise ValueError(
            "unsupported discrete-history transform; expected one of: lambda, kappa, delta, early-burst"
        )
    return resolved


def _validate_discrete_history_transform_request(
    *,
    transform: str | None,
    transform_parameter_value: float | None,
) -> None:
    if transform is None:
        if transform_parameter_value is not None:
            raise ValueError(
                "transform_parameter_value requires a discrete-history transform"
            )
        return
    if transform_parameter_value is None:
        raise ValueError(
            "transform_parameter_value is required when a discrete-history transform is supplied"
        )
    if transform == "lambda" and not 0.0 <= transform_parameter_value <= 1.0:
        raise ValueError("lambda transform_parameter_value must lie within [0, 1]")
    if transform == "kappa" and not 0.0 <= transform_parameter_value <= 1.0:
        raise ValueError("kappa transform_parameter_value must lie within [0, 1]")
    if transform == "delta" and transform_parameter_value <= 0.0:
        raise ValueError("delta transform_parameter_value must be positive")


def _discrete_history_transform_mode_name(transform: str) -> str:
    if transform == "lambda":
        return "pagel-lambda"
    if transform == "kappa":
        return "pagel-kappa"
    if transform == "delta":
        return "pagel-delta"
    if transform == "early-burst":
        return "early-burst"
    raise ValueError(
        "unsupported discrete-history transform; expected one of: lambda, kappa, delta, early-burst"
    )


def _discrete_history_transform_parameter_name(transform: str) -> str:
    if transform == "early-burst":
        return "a"
    return transform


def _transform_discrete_history_tree(
    tree: PhyloTree,
    *,
    transform: str,
    transform_parameter_value: float,
) -> PhyloTree:
    from bijux_phylogenetics.comparative.evolutionary_modes import (
        transform_tree_for_evolutionary_mode,
    )

    parameter_value = (
        -transform_parameter_value
        if transform == "early-burst"
        else transform_parameter_value
    )
    return transform_tree_for_evolutionary_mode(
        tree,
        mode=_discrete_history_transform_mode_name(transform),
        parameter_value=parameter_value,
        sigsq=1.0,
    )


def _summarize_discrete_history_collection(
    simulations: list[DiscreteTraitSimulationReport],
    *,
    states: tuple[str, ...],
):
    from .._statistics import _mean, _round_float
    from ..contracts import DiscreteHistorySummaryRow

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
    transform: str | None = None,
    transform_parameter_value: float | None = None,
    replicates: int = 1,
    seed: int = 1,
):
    from ..contracts import DiscreteHistorySimulationCollectionReport
    from .policy import (
        _normalize_discrete_states,
        _normalize_root_state_probabilities,
    )

    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    normalized_states = _normalize_discrete_states(states)
    resolved_transform = _resolve_discrete_history_transform_name(transform)
    _validate_discrete_history_transform_request(
        transform=resolved_transform,
        transform_parameter_value=transform_parameter_value,
    )
    normalized_rates = _normalize_rate_rows(
        states=normalized_states,
        rate_rows=rate_rows,
    )
    normalized_root_probabilities = _normalize_root_state_probabilities(
        states=normalized_states,
        root_state=root_state,
        root_state_probabilities=root_state_probabilities,
    )
    base_tree = load_tree(tree_path)
    working_tree = (
        base_tree
        if resolved_transform is None
        else _transform_discrete_history_tree(
            base_tree,
            transform=resolved_transform,
            transform_parameter_value=float(transform_parameter_value),
        )
    )
    transform_parameter_name = (
        None
        if resolved_transform is None
        else _discrete_history_transform_parameter_name(resolved_transform)
    )
    simulations = [
        _simulate_discrete_history_once(
            working_tree,
            model="rate-matrix-discrete-history",
            tree_path=tree_path,
            states=normalized_states,
            rate_rows=normalized_rates,
            transition_rate=None,
            fixed_root_state=root_state,
            root_state_probabilities=normalized_root_probabilities,
            transform_name=resolved_transform,
            transform_parameter_name=transform_parameter_name,
            transform_parameter_value=transform_parameter_value,
            seed=seed + index - 1,
        )
        for index in range(1, replicates + 1)
    ]
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
        tip_count=base_tree.tip_count,
        branch_count=sum(1 for _ in base_tree.iter_nodes()) - 1,
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
        transform_name=resolved_transform,
        transform_parameter_name=transform_parameter_name,
        transform_parameter_value=transform_parameter_value,
    )


def write_discrete_history_tip_truth_table(path: Path, report) -> Path:
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


def write_discrete_history_node_truth_table(path: Path, report) -> Path:
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


def write_discrete_history_branch_truth_table(path: Path, report) -> Path:
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


def write_discrete_history_event_table(path: Path, report) -> Path:
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


def write_discrete_history_segment_table(path: Path, report) -> Path:
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


def write_discrete_history_summary_table(path: Path, report) -> Path:
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
