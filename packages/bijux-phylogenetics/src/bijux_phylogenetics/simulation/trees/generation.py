from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import random

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.newick import (
    dumps_newick,
    write_newick,
    write_newick_tree_set,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.trees.inspection import summarize_tree_shape_from_tree

from .._statistics import (
    _mean,
    _median,
    _population_standard_deviation,
    _round_float,
)
from ..contracts import (
    CoalescentSkylineSummaryRow,
    CoalescentWaitingTimeSummaryRow,
    SimulatedTreeRecord,
    TreeSimulationEnvelopeMetric,
    TreeSimulationReport,
)

_DEFAULT_COALESCENT_WAITING_TIME_TOLERANCE = 0.2


@dataclass(slots=True)
class _Lineage:
    node: TreeNode
    start_time: float
    is_root: bool
    extant: bool = True


@dataclass(frozen=True, slots=True)
class _CoalescentWaitingInterval:
    lineage_count: int
    waiting_time: float


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
    return [
        float(node.branch_length or 0.0)
        for node in _iter_non_root_nodes_preorder(tree.root)
    ]


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
    coalescent_waiting_time_tolerance: float | None = None,
    coalescent_waiting_time_rows: list[CoalescentWaitingTimeSummaryRow] | None = None,
    coalescent_skyline_rows: list[CoalescentSkylineSummaryRow] | None = None,
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
        coalescent_waiting_time_tolerance=coalescent_waiting_time_tolerance,
        coalescent_waiting_time_rows=list(coalescent_waiting_time_rows or []),
        coalescent_skyline_rows=list(coalescent_skyline_rows or []),
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


def _simulate_random_tree_topology(
    node: TreeNode, tip_count: int, rng: random.Random
) -> None:
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


def _summarize_coalescent_waiting_times(
    waiting_intervals: list[_CoalescentWaitingInterval],
    *,
    population_size: float,
    tolerance: float,
) -> list[CoalescentWaitingTimeSummaryRow]:
    waiting_times_by_lineage_count: dict[int, list[float]] = {}
    for interval in waiting_intervals:
        waiting_times_by_lineage_count.setdefault(interval.lineage_count, []).append(
            interval.waiting_time
        )
    rows: list[CoalescentWaitingTimeSummaryRow] = []
    for lineage_count in sorted(waiting_times_by_lineage_count, reverse=True):
        waiting_times = waiting_times_by_lineage_count[lineage_count]
        coalescent_rate = lineage_count * (lineage_count - 1) / (2.0 * population_size)
        expected_waiting_time = 1.0 / coalescent_rate
        mean_waiting_time = _mean(waiting_times)
        absolute_error = abs(mean_waiting_time - expected_waiting_time)
        relative_error = (
            0.0
            if expected_waiting_time == 0.0
            else absolute_error / expected_waiting_time
        )
        rows.append(
            CoalescentWaitingTimeSummaryRow(
                lineage_count=lineage_count,
                coalescent_rate=_round_float(coalescent_rate),
                expected_waiting_time=_round_float(expected_waiting_time),
                observation_count=len(waiting_times),
                mean_waiting_time=mean_waiting_time,
                standard_deviation=_population_standard_deviation(waiting_times),
                minimum_waiting_time=_round_float(min(waiting_times)),
                median_waiting_time=_median(waiting_times),
                maximum_waiting_time=_round_float(max(waiting_times)),
                absolute_error=_round_float(absolute_error),
                relative_error=_round_float(relative_error),
                within_tolerance=relative_error <= tolerance,
            )
        )
    return rows


def _summarize_coalescent_skyline(
    waiting_time_rows: list[CoalescentWaitingTimeSummaryRow],
) -> list[CoalescentSkylineSummaryRow]:
    rows: list[CoalescentSkylineSummaryRow] = []
    for row in waiting_time_rows:
        effective_population_size_estimate = _round_float(
            row.mean_waiting_time * row.lineage_count * (row.lineage_count - 1) / 2.0
        )
        rows.append(
            CoalescentSkylineSummaryRow(
                interval=f"{row.lineage_count}->{row.lineage_count - 1}",
                lineage_count=row.lineage_count,
                duration=row.mean_waiting_time,
                effective_population_size_estimate=effective_population_size_estimate,
                observation_count=row.observation_count,
                relative_error=row.relative_error,
                uncertainty_flag=("low" if row.within_tolerance else "high"),
            )
        )
    return rows


def _simulate_coalescent_tree_once(
    *,
    tip_count: int,
    population_size: float,
    seed: int,
    taxon_prefix: str,
) -> tuple[PhyloTree, list[_CoalescentWaitingInterval]]:
    if population_size <= 0.0:
        raise ValueError(f"population_size must be positive, got {population_size}")

    rng = random.Random(seed)  # nosec B311
    node_heights = []
    cumulative_height = 0.0
    waiting_intervals: list[_CoalescentWaitingInterval] = []
    for lineage_count in range(tip_count, 1, -1):
        waiting_time = rng.expovariate(
            lineage_count * (lineage_count - 1) / (2.0 * population_size)
        )
        waiting_intervals.append(
            _CoalescentWaitingInterval(
                lineage_count=lineage_count,
                waiting_time=_round_float(waiting_time),
            )
        )
        cumulative_height += waiting_time
        node_heights.append(cumulative_height)
    pool: list[tuple[TreeNode, float]] = [(TreeNode(), 0.0) for _ in range(tip_count)]
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
    return tree, waiting_intervals


def simulate_coalescent_trees(
    *,
    tree_count: int,
    tip_count: int,
    population_size: float = 1.0,
    waiting_time_tolerance: float = _DEFAULT_COALESCENT_WAITING_TIME_TOLERANCE,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
) -> tuple[list[PhyloTree], TreeSimulationReport]:
    """Simulate rooted trees under a basic Kingman-style coalescent."""
    _validate_tree_count(tree_count, tip_count)
    if waiting_time_tolerance < 0.0:
        raise ValueError(
            f"waiting_time_tolerance must be nonnegative, got {waiting_time_tolerance}"
        )
    simulated = [
        _simulate_coalescent_tree_once(
            tip_count=tip_count,
            population_size=population_size,
            seed=seed + index - 1,
            taxon_prefix=taxon_prefix,
        )
        for index in range(1, tree_count + 1)
    ]
    trees = [tree for tree, _waiting_intervals in simulated]
    waiting_intervals = [
        interval
        for _tree, tree_waiting_intervals in simulated
        for interval in tree_waiting_intervals
    ]
    waiting_time_rows = _summarize_coalescent_waiting_times(
        waiting_intervals,
        population_size=population_size,
        tolerance=waiting_time_tolerance,
    )
    skyline_rows = _summarize_coalescent_skyline(waiting_time_rows)
    return trees, _build_tree_simulation_report(
        model="coalescent",
        tree_count=tree_count,
        tip_count=tip_count,
        seed=seed,
        trees=trees,
        branch_length_model="coalescent-waiting-times",
        population_size=population_size,
        coalescent_waiting_time_tolerance=waiting_time_tolerance,
        coalescent_waiting_time_rows=waiting_time_rows,
        coalescent_skyline_rows=skyline_rows,
    )


def simulate_coalescent_tree(
    *,
    tip_count: int,
    population_size: float = 1.0,
    waiting_time_tolerance: float = _DEFAULT_COALESCENT_WAITING_TIME_TOLERANCE,
    seed: int = 1,
    taxon_prefix: str = "Taxon",
) -> tuple[PhyloTree, TreeSimulationReport]:
    """Simulate one rooted coalescent tree with a structured report."""
    trees, report = simulate_coalescent_trees(
        tree_count=1,
        tip_count=tip_count,
        population_size=population_size,
        waiting_time_tolerance=waiting_time_tolerance,
        seed=seed,
        taxon_prefix=taxon_prefix,
    )
    return trees[0], report


def write_tree_set(path: Path, trees: list[PhyloTree]) -> Path:
    """Write a list of simulated trees as one canonical Newick tree per line."""
    return write_newick_tree_set(path, trees)


def write_simulated_tree(path: Path, tree: PhyloTree) -> Path:
    """Write one simulated tree as canonical Newick."""
    return write_newick(path, tree)


def write_tree_simulation_record_table(
    path: Path, report: TreeSimulationReport
) -> Path:
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


def write_coalescent_waiting_time_table(
    path: Path,
    report: TreeSimulationReport,
) -> Path:
    """Write one row per lineage-count waiting-time summary from a coalescent report."""
    return write_taxon_rows(
        path,
        columns=[
            "lineage_count",
            "coalescent_rate",
            "expected_waiting_time",
            "observation_count",
            "mean_waiting_time",
            "standard_deviation",
            "minimum_waiting_time",
            "median_waiting_time",
            "maximum_waiting_time",
            "absolute_error",
            "relative_error",
            "within_tolerance",
            "waiting_time_tolerance",
        ],
        rows=[
            {
                "lineage_count": str(row.lineage_count),
                "coalescent_rate": format(row.coalescent_rate, ".15g"),
                "expected_waiting_time": format(row.expected_waiting_time, ".15g"),
                "observation_count": str(row.observation_count),
                "mean_waiting_time": format(row.mean_waiting_time, ".15g"),
                "standard_deviation": format(row.standard_deviation, ".15g"),
                "minimum_waiting_time": format(row.minimum_waiting_time, ".15g"),
                "median_waiting_time": format(row.median_waiting_time, ".15g"),
                "maximum_waiting_time": format(row.maximum_waiting_time, ".15g"),
                "absolute_error": format(row.absolute_error, ".15g"),
                "relative_error": format(row.relative_error, ".15g"),
                "within_tolerance": str(row.within_tolerance).lower(),
                "waiting_time_tolerance": (
                    ""
                    if report.coalescent_waiting_time_tolerance is None
                    else format(report.coalescent_waiting_time_tolerance, ".15g")
                ),
            }
            for row in report.coalescent_waiting_time_rows
        ],
    )


def write_coalescent_skyline_table(
    path: Path,
    report: TreeSimulationReport,
) -> Path:
    """Write one skyline interval row per lineage-count transition from a coalescent report."""
    return write_taxon_rows(
        path,
        columns=[
            "interval",
            "lineage_count",
            "duration",
            "effective_population_size_estimate",
            "observation_count",
            "relative_error",
            "uncertainty_flag",
        ],
        rows=[
            {
                "interval": row.interval,
                "lineage_count": str(row.lineage_count),
                "duration": format(row.duration, ".15g"),
                "effective_population_size_estimate": format(
                    row.effective_population_size_estimate,
                    ".15g",
                ),
                "observation_count": str(row.observation_count),
                "relative_error": format(row.relative_error, ".15g"),
                "uncertainty_flag": row.uncertainty_flag,
            }
            for row in report.coalescent_skyline_rows
        ],
    )
