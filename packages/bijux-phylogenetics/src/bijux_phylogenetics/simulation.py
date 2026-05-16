from __future__ import annotations

from dataclasses import dataclass, field
from math import exp, sqrt
from pathlib import Path
import random
from statistics import median

from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.core.alignment import AlignmentRecord
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.io.fasta import write_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick, write_newick, write_newick_tree_set
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.tree_shape import summarize_tree_shape_from_tree


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
    alpha: float | None
    theta: float | None
    rate_change: float | None
    traits: list[SimulatedContinuousTrait]
    node_values: list[SimulatedContinuousNode]


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


@dataclass(slots=True)
class _Lineage:
    node: TreeNode
    start_time: float
    is_root: bool
    extant: bool = True


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


def _validate_tree_count(tree_count: int, tip_count: int) -> None:
    if tree_count < 1:
        raise ValueError(f"tree_count must be at least 1, got {tree_count}")
    if tip_count < 2:
        raise ValueError(f"tip_count must be at least 2, got {tip_count}")


def _finalize_branch(lineage: _Lineage, event_time: float) -> None:
    if not lineage.is_root:
        lineage.node.branch_length = round(event_time - lineage.start_time, 15)


def _prune_extinct_subtree(node: TreeNode, *, is_root: bool = False) -> TreeNode | None:
    if node.is_leaf():
        return node if node.name is not None else None

    kept_children = [
        pruned_child
        for child in node.children
        if (pruned_child := _prune_extinct_subtree(child)) is not None
    ]
    if not kept_children:
        return None
    if len(kept_children) == 1:
        promoted = kept_children[0]
        if not is_root and node.branch_length is not None:
            promoted.branch_length = round(
                (promoted.branch_length or 0.0) + node.branch_length, 15
            )
        return promoted
    node.replace_children(kept_children)
    return node


def _finalize_extant_leaves(lineages: list[_Lineage], *, final_time: float) -> None:
    for index, lineage in enumerate(lineages, start=1):
        _finalize_branch(lineage, final_time)
        lineage.node.name = f"__extant_{index}"


def _label_tree_leaves(tree: PhyloTree, *, taxon_prefix: str) -> None:
    for index, leaf in enumerate(tree.iter_leaves(), start=1):
        leaf.name = f"{taxon_prefix}{index}"


def _label_tree_leaves_randomized(
    tree: PhyloTree,
    *,
    taxon_prefix: str,
    rng: random.Random,
) -> None:
    labels = [f"{taxon_prefix}{index}" for index in range(1, tree.tip_count + 1)]
    rng.shuffle(labels)
    for leaf, label in zip(tree.iter_leaves(), labels, strict=True):
        leaf.name = label


def _iter_non_root_nodes_preorder(node: TreeNode):
    for child in node.children:
        yield child
        yield from _iter_non_root_nodes_preorder(child)


def _branch_lengths(tree: PhyloTree) -> list[float]:
    return [float(node.branch_length or 0.0) for node in _iter_non_root_nodes_preorder(tree.root)]


def _simulation_envelope_metric(
    metric: str,
    sample_scope: str,
    values: list[float],
) -> TreeSimulationEnvelopeMetric:
    return TreeSimulationEnvelopeMetric(
        metric=metric,
        sample_scope=sample_scope,
        observation_count=len(values),
        mean=_mean(values),
        standard_deviation=_population_standard_deviation(values),
        minimum=_round_float(min(values)),
        median=_median(values),
        maximum=_round_float(max(values)),
    )


def _build_simulated_tree_record(
    tree: PhyloTree,
    *,
    index: int,
) -> tuple[SimulatedTreeRecord, list[float]]:
    branch_lengths = _branch_lengths(tree)
    shape = summarize_tree_shape_from_tree(
        tree,
        source_path=Path("simulated-tree.nwk"),
        tree_index=index,
    )
    total_branch_length = _round_float(sum(branch_lengths))
    normalized_colless = (
        0.0
        if shape.normalized_colless_imbalance is None
        else _round_float(shape.normalized_colless_imbalance)
    )
    record = SimulatedTreeRecord(
        index=index,
        newick=dumps_newick(tree),
        tree_height_branch_length=_round_float(shape.tree_height_branch_length or 0.0),
        total_branch_length=total_branch_length,
        mean_branch_length=_mean(branch_lengths),
        median_branch_length=_median(branch_lengths),
        minimum_branch_length=_round_float(min(branch_lengths)),
        maximum_branch_length=_round_float(max(branch_lengths)),
        cherry_count=shape.cherry_count,
        sackin_imbalance_index=shape.sackin_imbalance_index,
        normalized_colless_imbalance=normalized_colless,
    )
    return record, branch_lengths


def _build_tree_simulation_report(
    *,
    model: str,
    tree_count: int,
    tip_count: int,
    seed: int,
    trees: list[PhyloTree],
    branch_length_model: str | None,
    birth_rate: float | None = None,
    death_rate: float | None = None,
    population_size: float | None = None,
) -> TreeSimulationReport:
    records: list[SimulatedTreeRecord] = []
    pooled_branch_lengths: list[float] = []
    for index, tree in enumerate(trees, start=1):
        record, branch_lengths = _build_simulated_tree_record(tree, index=index)
        records.append(record)
        pooled_branch_lengths.extend(branch_lengths)
    envelope_metrics = [
        _simulation_envelope_metric(
            "tree_height_branch_length",
            "tree",
            [record.tree_height_branch_length for record in records],
        ),
        _simulation_envelope_metric(
            "total_branch_length",
            "tree",
            [record.total_branch_length for record in records],
        ),
        _simulation_envelope_metric(
            "branch_length",
            "edge",
            pooled_branch_lengths,
        ),
        _simulation_envelope_metric(
            "cherry_count",
            "tree",
            [float(record.cherry_count) for record in records],
        ),
        _simulation_envelope_metric(
            "sackin_imbalance_index",
            "tree",
            [float(record.sackin_imbalance_index) for record in records],
        ),
        _simulation_envelope_metric(
            "normalized_colless_imbalance",
            "tree",
            [record.normalized_colless_imbalance for record in records],
        ),
    ]
    return TreeSimulationReport(
        model=model,
        tree_count=tree_count,
        tip_count=tip_count,
        seed=seed,
        records=records,
        branch_length_model=branch_length_model,
        birth_rate=birth_rate,
        death_rate=death_rate,
        population_size=population_size,
        rooted=all(tree.rooted for tree in trees),
        binary=all(tree.internal_node_count == tree.tip_count - 1 for tree in trees),
        pooled_branch_count=len(pooled_branch_lengths),
        envelope_metrics=envelope_metrics,
    )


def _simulate_birth_death_tree_once(
    *,
    tip_count: int,
    birth_rate: float,
    death_rate: float,
    seed: int,
    taxon_prefix: str,
) -> PhyloTree:
    if birth_rate <= 0.0:
        raise ValueError(f"birth_rate must be positive, got {birth_rate}")
    if death_rate < 0.0:
        raise ValueError(f"death_rate must be nonnegative, got {death_rate}")

    rng = random.Random(seed)  # nosec B311
    for attempt in range(1, 129):
        root = TreeNode()
        extant = [_Lineage(node=root, start_time=0.0, is_root=True)]
        absolute_time = 0.0
        while 0 < len(extant) < tip_count:
            total_rate = len(extant) * (birth_rate + death_rate)
            absolute_time += rng.expovariate(total_rate)
            selected_index = rng.randrange(len(extant))
            lineage = extant.pop(selected_index)
            event_is_birth = len(extant) == 0 or rng.random() < (
                birth_rate / (birth_rate + death_rate)
            )
            _finalize_branch(lineage, absolute_time)
            if event_is_birth:
                left = TreeNode()
                right = TreeNode()
                lineage.node.replace_children([left, right])
                extant.extend(
                    [
                        _Lineage(node=left, start_time=absolute_time, is_root=False),
                        _Lineage(node=right, start_time=absolute_time, is_root=False),
                    ]
                )
        if len(extant) != tip_count:
            rng.seed(seed + attempt)
            continue
        _finalize_extant_leaves(extant, final_time=absolute_time)
        pruned_root = _prune_extinct_subtree(root, is_root=True)
        if pruned_root is None:
            rng.seed(seed + attempt)
            continue
        pruned_root.branch_length = None
        tree = PhyloTree(root=pruned_root, source_format="newick")
        tree.rooted = True
        _label_tree_leaves(tree, taxon_prefix=taxon_prefix)
        tree.refresh()
        return tree
    raise ValueError(
        "birth-death simulation failed to retain the requested number of extant tips after 128 attempts"
    )


def _simulate_random_tree_topology(node: TreeNode, tip_count: int, rng: random.Random) -> None:
    if tip_count == 1:
        return
    if tip_count == 2:
        node.replace_children([TreeNode(), TreeNode()])
        return
    left_tip_count = rng.randrange(1, tip_count)
    right_tip_count = tip_count - left_tip_count
    left = TreeNode()
    right = TreeNode()
    node.replace_children([left, right])
    _simulate_random_tree_topology(left, left_tip_count, rng)
    _simulate_random_tree_topology(right, right_tip_count, rng)


def _simulate_random_tree_once(
    *,
    tip_count: int,
    seed: int,
    taxon_prefix: str,
    branch_length_model: str,
) -> PhyloTree:
    if branch_length_model != "uniform":
        raise ValueError(
            "random-tree simulation currently supports only the 'uniform' branch-length model"
        )
    rng = random.Random(seed)  # nosec B311
    root = TreeNode()
    _simulate_random_tree_topology(root, tip_count, rng)
    tree = PhyloTree(root=root, source_format="newick")
    tree.rooted = True
    _label_tree_leaves_randomized(tree, taxon_prefix=taxon_prefix, rng=rng)
    tree.refresh()
    for node in _iter_non_root_nodes_preorder(tree.root):
        node.branch_length = _round_float(rng.random())
    return tree


def simulate_birth_death_trees(
    *,
    tree_count: int,
    tip_count: int,
    birth_rate: float = 1.0,
    death_rate: float = 0.25,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
) -> tuple[list[PhyloTree], TreeSimulationReport]:
    """Simulate rooted extant trees under a simple birth-death process."""
    _validate_tree_count(tree_count, tip_count)
    trees = [
        _simulate_birth_death_tree_once(
            tip_count=tip_count,
            birth_rate=birth_rate,
            death_rate=death_rate,
            seed=seed + index - 1,
            taxon_prefix=taxon_prefix,
        )
        for index in range(1, tree_count + 1)
    ]
    return trees, _build_tree_simulation_report(
        model="birth-death",
        tree_count=tree_count,
        tip_count=tip_count,
        seed=seed,
        trees=trees,
        branch_length_model="birth-death",
        birth_rate=birth_rate,
        death_rate=death_rate,
    )


def simulate_random_trees(
    *,
    tree_count: int,
    tip_count: int,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
    branch_length_model: str = "uniform",
) -> tuple[list[PhyloTree], TreeSimulationReport]:
    """Simulate rooted binary random trees with ape-style default branch lengths."""
    _validate_tree_count(tree_count, tip_count)
    trees = [
        _simulate_random_tree_once(
            tip_count=tip_count,
            seed=seed + index - 1,
            taxon_prefix=taxon_prefix,
            branch_length_model=branch_length_model,
        )
        for index in range(1, tree_count + 1)
    ]
    return trees, _build_tree_simulation_report(
        model="random-tree",
        tree_count=tree_count,
        tip_count=tip_count,
        seed=seed,
        trees=trees,
        branch_length_model=branch_length_model,
    )


def simulate_random_tree(
    *,
    tip_count: int,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
    branch_length_model: str = "uniform",
) -> tuple[PhyloTree, TreeSimulationReport]:
    """Simulate one rooted binary random tree with a structured report."""
    trees, report = simulate_random_trees(
        tree_count=1,
        tip_count=tip_count,
        seed=seed,
        taxon_prefix=taxon_prefix,
        branch_length_model=branch_length_model,
    )
    return trees[0], report


def _choose_two_indices(rng: random.Random, count: int) -> tuple[int, int]:
    left = rng.randrange(count)
    right = rng.randrange(count - 1)
    if right >= left:
        right += 1
    return tuple(sorted((left, right)))


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
    return {
        node.name: (
            round(float(node_values[node_signature(node)]), 15)
            if isinstance(node_values[node_signature(node)], float)
            else node_values[node_signature(node)]
        )
        for node in tree.iter_leaves()
        if node.name is not None
    }


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
            raise ValueError(f"unknown source_state '{row.source_state}' in rate matrix")
        if row.target_state not in state_set:
            raise ValueError(f"unknown target_state '{row.target_state}' in rate matrix")
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
    return {state: _round_float(value / total) for state, value in probabilities.items()}


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
        source_state: {target_state: 0.0 for target_state in states if target_state != source_state}
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
        branch_transition_counts = {label: 0.0 for label in transition_labels}
        branch_state_totals = {state: 0.0 for state in states}
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


def _simulate_coalescent_tree_once(
    *,
    tip_count: int,
    population_size: float,
    seed: int,
    taxon_prefix: str,
) -> PhyloTree:
    if population_size <= 0.0:
        raise ValueError(f"population_size must be positive, got {population_size}")

    rng = random.Random(seed)  # nosec B311
    node_heights = []
    cumulative_height = 0.0
    for lineage_count in range(tip_count, 1, -1):
        cumulative_height += rng.expovariate(
            lineage_count * (lineage_count - 1) / (2.0 * population_size)
        )
        node_heights.append(cumulative_height)
    pool: list[tuple[TreeNode, float]] = [
        (TreeNode(), 0.0) for _ in range(tip_count)
    ]
    for node_height in node_heights:
        left_index, right_index = _choose_two_indices(rng, len(pool))
        right_node, right_height = pool.pop(right_index)
        left_node, left_height = pool.pop(left_index)
        left_node.branch_length = _round_float(node_height - left_height)
        right_node.branch_length = _round_float(node_height - right_height)
        pool.append((TreeNode(children=[left_node, right_node]), node_height))
    root, _root_height = pool[0]
    root.branch_length = None
    tree = PhyloTree(root=root, source_format="newick")
    tree.rooted = True
    _label_tree_leaves_randomized(tree, taxon_prefix=taxon_prefix, rng=rng)
    return tree


def simulate_coalescent_trees(
    *,
    tree_count: int,
    tip_count: int,
    population_size: float = 1.0,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
) -> tuple[list[PhyloTree], TreeSimulationReport]:
    """Simulate rooted trees under a basic Kingman-style coalescent."""
    _validate_tree_count(tree_count, tip_count)
    trees = [
        _simulate_coalescent_tree_once(
            tip_count=tip_count,
            population_size=population_size,
            seed=seed + index - 1,
            taxon_prefix=taxon_prefix,
        )
        for index in range(1, tree_count + 1)
    ]
    return trees, _build_tree_simulation_report(
        model="coalescent",
        tree_count=tree_count,
        tip_count=tip_count,
        seed=seed,
        trees=trees,
        branch_length_model="coalescent-waiting-times",
        population_size=population_size,
    )


def simulate_coalescent_tree(
    *,
    tip_count: int,
    population_size: float = 1.0,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
) -> tuple[PhyloTree, TreeSimulationReport]:
    """Simulate one rooted coalescent tree with a structured report."""
    trees, report = simulate_coalescent_trees(
        tree_count=1,
        tip_count=tip_count,
        population_size=population_size,
        seed=seed,
        taxon_prefix=taxon_prefix,
    )
    return trees[0], report


def simulate_brownian_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float = 1.0,
    seed: int = 1,
) -> ContinuousTraitSimulationReport:
    """Simulate one continuous tip trait under Brownian motion."""
    if sigma < 0.0:
        raise ValueError(f"sigma must be nonnegative, got {sigma}")
    tree = load_tree(tree_path)
    rng = random.Random(seed)  # nosec B311
    node_values = _iter_node_trait_values(
        tree,
        root_state=root_state,
        propagate=lambda state, branch_length: (
            state + rng.gauss(0.0, sigma * sqrt(branch_length))
        ),
    )
    return _build_continuous_trait_simulation_report(
        tree=tree,
        tree_path=tree_path,
        model="brownian-motion",
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        alpha=None,
        theta=None,
        rate_change=None,
        node_values=node_values,
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
    """Simulate one continuous tip trait under an OU process."""
    if sigma < 0.0:
        raise ValueError(f"sigma must be nonnegative, got {sigma}")
    if alpha < 0.0:
        raise ValueError(f"alpha must be nonnegative, got {alpha}")
    tree = load_tree(tree_path)
    rng = random.Random(seed)  # nosec B311

    def propagate(state: float, branch_length: float) -> float:
        if branch_length == 0.0:
            return state
        if alpha == 0.0:
            return state + rng.gauss(0.0, sigma * sqrt(branch_length))
        mean = theta + (state - theta) * exp(-alpha * branch_length)
        variance = (
            (sigma**2) * (1.0 - exp(-2.0 * alpha * branch_length)) / (2.0 * alpha)
        )
        return mean + rng.gauss(0.0, sqrt(max(variance, 0.0)))

    node_values = _iter_node_trait_values(
        tree, root_state=root_state, propagate=propagate
    )
    return _build_continuous_trait_simulation_report(
        tree=tree,
        tree_path=tree_path,
        model="ornstein-uhlenbeck",
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        alpha=alpha,
        theta=theta,
        rate_change=None,
        node_values=node_values,
    )


def simulate_early_burst_traits(
    tree_path: Path,
    *,
    root_state: float = 0.0,
    sigma: float = 1.0,
    rate_change: float = 1.0,
    seed: int = 1,
) -> ContinuousTraitSimulationReport:
    """Simulate one continuous tip trait under an early-burst branch-rate process."""
    if sigma < 0.0:
        raise ValueError(f"sigma must be nonnegative, got {sigma}")
    if rate_change < 0.0:
        raise ValueError(f"rate_change must be nonnegative, got {rate_change}")
    tree = load_tree(tree_path)
    from bijux_phylogenetics.comparative.evolutionary_modes import (
        transform_tree_for_evolutionary_mode,
    )

    transformed_tree = transform_tree_for_evolutionary_mode(
        tree,
        mode="early-burst",
        parameter_value=rate_change,
    )
    rng = random.Random(seed)  # nosec B311
    node_values = _iter_node_trait_values(
        transformed_tree,
        root_state=root_state,
        propagate=lambda state, branch_length: (
            state + rng.gauss(0.0, sigma * sqrt(branch_length))
        ),
    )
    return _build_continuous_trait_simulation_report(
        tree=transformed_tree,
        tree_path=tree_path,
        model="early-burst",
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        alpha=None,
        theta=None,
        rate_change=rate_change,
        node_values=node_values,
    )


def _build_continuous_trait_simulation_report(
    *,
    tree: PhyloTree,
    tree_path: Path,
    model: str,
    seed: int,
    root_state: float,
    sigma: float,
    alpha: float | None,
    theta: float | None,
    rate_change: float | None,
    node_values: dict[str, float],
) -> ContinuousTraitSimulationReport:
    values = _tip_values_from_node_map(tree, node_values)
    return ContinuousTraitSimulationReport(
        model=model,
        tree_path=tree_path,
        tip_count=tree.tip_count,
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        alpha=alpha,
        theta=theta,
        rate_change=rate_change,
        traits=[
            SimulatedContinuousTrait(taxon=taxon, value=value)
            for taxon, value in sorted(values.items())
        ],
        node_values=[
            SimulatedContinuousNode(
                node=node_signature(node),
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=node_descendant_taxa(node),
                value=float(format(node_values[node_signature(node)], ".15g")),
            )
            for node in tree.iter_nodes()
        ],
    )


def simulate_discrete_traits(
    tree_path: Path,
    *,
    states: list[str],
    transition_rate: float = 1.0,
    root_state: str | None = None,
    seed: int = 1,
) -> DiscreteTraitSimulationReport:
    """Simulate one discrete tip trait over a tree using symmetric jump changes."""
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

    def visit(node: TreeNode, state: str) -> None:
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


def _simulate_alignment_records(
    tree_path: Path,
    *,
    alphabet: tuple[str, ...],
    model: str,
    sequence_length: int,
    substitution_rate: float,
    seed: int,
) -> AlignmentSimulationReport:
    if sequence_length < 1:
        raise ValueError(f"sequence_length must be at least 1, got {sequence_length}")
    if substitution_rate < 0.0:
        raise ValueError(
            f"substitution_rate must be nonnegative, got {substitution_rate}"
        )
    tree = load_tree(tree_path)
    rng = random.Random(seed)  # nosec B311
    root_sequence = "".join(rng.choice(alphabet) for _ in range(sequence_length))

    def mutate_sequence(sequence: str, branch_length: float) -> str:
        residues: list[str] = []
        for residue in sequence:
            next_residue = residue
            for _ in range(_poisson_count(substitution_rate * branch_length, rng)):
                alternatives = [
                    candidate for candidate in alphabet if candidate != next_residue
                ]
                next_residue = rng.choice(alternatives)
            residues.append(next_residue)
        return "".join(residues)

    values = _iter_tip_trait_values(
        tree,
        root_state=root_sequence,
        propagate=lambda state, branch_length: mutate_sequence(state, branch_length),
    )
    return AlignmentSimulationReport(
        model=model,
        tree_path=tree_path,
        tip_count=tree.tip_count,
        seed=seed,
        sequence_length=sequence_length,
        substitution_rate=substitution_rate,
        inferred_alphabet="dna" if alphabet == ("A", "C", "G", "T") else "protein",
        records=[
            AlignmentRecord(identifier=taxon, sequence=sequence)
            for taxon, sequence in sorted(values.items())
        ],
    )


def simulate_dna_alignment(
    tree_path: Path,
    *,
    sequence_length: int,
    substitution_rate: float = 1.0,
    seed: int = 1,
) -> AlignmentSimulationReport:
    """Simulate a DNA alignment along a rooted tree under a simple JC-like process."""
    return _simulate_alignment_records(
        tree_path,
        alphabet=("A", "C", "G", "T"),
        model="jukes-cantor-like",
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
    """Simulate a protein alignment along a rooted tree under a symmetric exchange model."""
    return _simulate_alignment_records(
        tree_path,
        alphabet=tuple("ACDEFGHIKLMNPQRSTVWY"),
        model="symmetric-protein",
        sequence_length=sequence_length,
        substitution_rate=substitution_rate,
        seed=seed,
    )


def write_tree_set(path: Path, trees: list[PhyloTree]) -> Path:
    """Write a list of simulated trees as one canonical Newick tree per line."""
    return write_newick_tree_set(path, trees)


def write_simulated_tree(path: Path, tree: PhyloTree) -> Path:
    """Write one simulated tree as canonical Newick."""
    return write_newick(path, tree)


def write_tree_simulation_record_table(path: Path, report: TreeSimulationReport) -> Path:
    """Write one per-tree simulation metrics row for a governed tree simulation report."""
    return write_taxon_rows(
        path,
        columns=[
            "index",
            "tree_height_branch_length",
            "total_branch_length",
            "mean_branch_length",
            "median_branch_length",
            "minimum_branch_length",
            "maximum_branch_length",
            "cherry_count",
            "sackin_imbalance_index",
            "normalized_colless_imbalance",
            "newick",
        ],
        rows=[
            {
                "index": str(record.index),
                "tree_height_branch_length": format(
                    record.tree_height_branch_length, ".15g"
                ),
                "total_branch_length": format(record.total_branch_length, ".15g"),
                "mean_branch_length": format(record.mean_branch_length, ".15g"),
                "median_branch_length": format(record.median_branch_length, ".15g"),
                "minimum_branch_length": format(record.minimum_branch_length, ".15g"),
                "maximum_branch_length": format(record.maximum_branch_length, ".15g"),
                "cherry_count": str(record.cherry_count),
                "sackin_imbalance_index": str(record.sackin_imbalance_index),
                "normalized_colless_imbalance": format(
                    record.normalized_colless_imbalance, ".15g"
                ),
                "newick": record.newick,
            }
            for record in report.records
        ],
    )


def write_tree_simulation_envelope_table(
    path: Path,
    report: TreeSimulationReport,
) -> Path:
    """Write the aggregate simulation-envelope metrics for a governed tree simulation report."""
    return write_taxon_rows(
        path,
        columns=[
            "metric",
            "sample_scope",
            "observation_count",
            "mean",
            "standard_deviation",
            "minimum",
            "median",
            "maximum",
        ],
        rows=[
            {
                "metric": row.metric,
                "sample_scope": row.sample_scope,
                "observation_count": str(row.observation_count),
                "mean": format(row.mean, ".15g"),
                "standard_deviation": format(row.standard_deviation, ".15g"),
                "minimum": format(row.minimum, ".15g"),
                "median": format(row.median, ".15g"),
                "maximum": format(row.maximum, ".15g"),
            }
            for row in report.envelope_metrics
        ],
    )


def write_continuous_trait_table(
    path: Path, report: ContinuousTraitSimulationReport
) -> Path:
    """Write simulated continuous trait values as a taxon-keyed table."""
    return write_taxon_rows(
        path,
        columns=["taxon", "value"],
        rows=[
            {"taxon": row.taxon, "value": format(row.value, ".15g")}
            for row in report.traits
        ],
    )


def write_discrete_trait_table(
    path: Path, report: DiscreteTraitSimulationReport
) -> Path:
    """Write simulated discrete trait states as a taxon-keyed table."""
    return write_taxon_rows(
        path,
        columns=["taxon", "state"],
        rows=[{"taxon": row.taxon, "state": row.state} for row in report.traits],
    )


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
    """Write a simulated alignment as FASTA."""
    return write_fasta_alignment(path, report.records)
