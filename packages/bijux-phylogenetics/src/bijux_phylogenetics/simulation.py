from __future__ import annotations

from dataclasses import dataclass
from math import exp, sqrt
from pathlib import Path
import random

from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.core.alignment import AlignmentRecord
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.io.fasta import write_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.io.trees import load_tree


@dataclass(frozen=True, slots=True)
class SimulatedTreeRecord:
    index: int
    newick: str


@dataclass(slots=True)
class TreeSimulationReport:
    model: str
    tree_count: int
    tip_count: int
    seed: int
    records: list[SimulatedTreeRecord]


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
        _label_tree_leaves(tree, taxon_prefix=taxon_prefix)
        return tree
    raise ValueError(
        "birth-death simulation failed to retain the requested number of extant tips after 128 attempts"
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
    return trees, TreeSimulationReport(
        model="birth-death",
        tree_count=tree_count,
        tip_count=tip_count,
        seed=seed,
        records=[
            SimulatedTreeRecord(index=index, newick=dumps_newick(tree))
            for index, tree in enumerate(trees, start=1)
        ],
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
) -> str:
    if rate < 0.0:
        raise ValueError(f"rate must be nonnegative, got {rate}")
    next_state = state
    for _ in range(_poisson_count(rate * branch_length, rng)):
        alternatives = [candidate for candidate in states if candidate != next_state]
        next_state = rng.choice(alternatives)
    return next_state


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
    return PhyloTree(root=root, source_format="newick")


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
    return trees, TreeSimulationReport(
        model="coalescent",
        tree_count=tree_count,
        tip_count=tip_count,
        seed=seed,
        records=[
            SimulatedTreeRecord(index=index, newick=dumps_newick(tree))
            for index, tree in enumerate(trees, start=1)
        ],
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
    values = _tip_values_from_node_map(tree, node_values)
    return ContinuousTraitSimulationReport(
        model="brownian-motion",
        tree_path=tree_path,
        tip_count=tree.tip_count,
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        alpha=None,
        theta=None,
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
    values = _tip_values_from_node_map(tree, node_values)
    return ContinuousTraitSimulationReport(
        model="ornstein-uhlenbeck",
        tree_path=tree_path,
        tip_count=tree.tip_count,
        seed=seed,
        root_state=root_state,
        sigma=sigma,
        alpha=alpha,
        theta=theta,
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
    node_values = _iter_node_trait_values(
        tree,
        root_state=starting_state,
        propagate=lambda state, branch_length: _simulate_symmetric_state_trajectory(
            state,
            branch_length=branch_length,
            rate=transition_rate,
            states=unique_states,
            rng=rng,
        ),
    )
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
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(f"{dumps_newick(tree)}\n" for tree in trees), encoding="utf-8"
    )
    return path


def write_simulated_tree(path: Path, tree: PhyloTree) -> Path:
    """Write one simulated tree as canonical Newick."""
    return write_newick(path, tree)


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
