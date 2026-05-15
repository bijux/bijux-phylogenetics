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


@dataclass(slots=True)
class DiscreteTraitSimulationReport:
    model: str
    tree_path: Path
    tip_count: int
    seed: int
    states: list[str]
    transition_rate: float
    root_state: str
    traits: list[SimulatedDiscreteTrait]
    node_states: list[SimulatedDiscreteNode]
    branch_histories: list[SimulatedDiscreteBranchHistory]


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
    node.children = kept_children
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
                lineage.node.children = [left, right]
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
        return tree
    raise ValueError(
        "birth-death simulation failed to retain the requested number of extant tips after 128 attempts"
    )


def _simulate_random_tree_topology(node: TreeNode, tip_count: int, rng: random.Random) -> None:
    if tip_count == 1:
        return
    if tip_count == 2:
        node.children = [TreeNode(), TreeNode()]
        return
    left_tip_count = rng.randrange(1, tip_count)
    right_tip_count = tip_count - left_tip_count
    left = TreeNode()
    right = TreeNode()
    node.children = [left, right]
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
) -> tuple[str, list[tuple[str, str]]]:
    if rate < 0.0:
        raise ValueError(f"rate must be nonnegative, got {rate}")
    next_state = state
    transitions: list[tuple[str, str]] = []
    for _ in range(_poisson_count(rate * branch_length, rng)):
        alternatives = [candidate for candidate in states if candidate != next_state]
        candidate = rng.choice(alternatives)
        transitions.append((next_state, candidate))
        next_state = candidate
    return next_state, transitions


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
    lineages = [
        _Lineage(
            node=TreeNode(name=f"{taxon_prefix}{index}"),
            start_time=0.0,
            is_root=False,
        )
        for index in range(1, tip_count + 1)
    ]
    absolute_time = 0.0
    while len(lineages) > 1:
        lineage_count = len(lineages)
        rate = lineage_count * (lineage_count - 1) / (2.0 * population_size)
        absolute_time += rng.expovariate(rate)
        left_index, right_index = _choose_two_indices(rng, len(lineages))
        right = lineages.pop(right_index)
        left = lineages.pop(left_index)
        _finalize_branch(left, absolute_time)
        _finalize_branch(right, absolute_time)
        parent = TreeNode(children=[left.node, right.node])
        lineages.append(_Lineage(node=parent, start_time=absolute_time, is_root=False))
    root = lineages[0].node
    root.branch_length = None
    tree = PhyloTree(root=root, source_format="newick")
    tree.rooted = True
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
    unique_states = tuple(dict.fromkeys(state for state in states if state))
    if len(unique_states) < 2:
        raise ValueError("states must contain at least two distinct non-empty states")
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
            child_state, transitions = _simulate_symmetric_state_trajectory(
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
                    changed=child_state != state,
                    event_count=len(transitions),
                    events=[
                        SimulatedDiscreteTransitionEvent(
                            parent_node=parent_signature,
                            child_node=child_signature,
                            source_state=source_state,
                            target_state=target_state,
                            event_index=index,
                        )
                        for index, (source_state, target_state) in enumerate(
                            transitions, start=1
                        )
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


def write_simulated_alignment(path: Path, report: AlignmentSimulationReport) -> Path:
    """Write a simulated alignment as FASTA."""
    return write_fasta_alignment(path, report.records)
