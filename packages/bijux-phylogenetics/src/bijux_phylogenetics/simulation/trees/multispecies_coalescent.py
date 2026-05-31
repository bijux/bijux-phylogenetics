from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import random

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
    assess_tree_ultrametricity,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .._statistics import _round_float
from ..contracts import (
    MultispeciesCoalescentBranchRow,
    MultispeciesCoalescentEventRow,
    MultispeciesCoalescentReport,
    MultispeciesCoalescentSampleRow,
)


@dataclass(frozen=True, slots=True)
class _ActiveGeneLineage:
    node: TreeNode
    age: float


@dataclass(frozen=True, slots=True)
class _SpeciesPopulation:
    species_branch: str
    branch_role: str
    descendant_species: list[str]
    population_size: float
    branch_start_age: float
    branch_end_age: float | None

    @property
    def branch_duration(self) -> float | None:
        if self.branch_end_age is None:
            return None
        return _round_float(self.branch_end_age - self.branch_start_age)

    @property
    def included_in_deep_coalescence_total(self) -> bool:
        return self.branch_role != "root-population"


def simulate_multispecies_coalescent_gene_tree(
    species_tree_path: Path,
    *,
    default_population_size: float = 1.0,
    sample_count_table_path: Path | None = None,
    population_size_table_path: Path | None = None,
    seed: int = 1,
) -> tuple[PhyloTree, MultispeciesCoalescentReport]:
    """Simulate one gene tree under a rooted ultrametric species tree."""
    if default_population_size <= 0.0:
        raise ValueError(
            f"default_population_size must be positive, got {default_population_size}"
        )
    species_tree = load_tree(species_tree_path)
    _validate_species_tree(species_tree, species_tree_path)
    ultrametric_report = assess_tree_ultrametricity(
        species_tree_path,
        tolerance=APE_ULTRAMETRIC_TOLERANCE,
    )
    if not ultrametric_report.ultrametric:
        raise ValueError(
            "multispecies coalescent simulation requires an ultrametric species tree"
        )
    node_age_by_id = _node_age_by_id(species_tree, root_age=ultrametric_report.root_age)
    sample_count_by_species = _load_sample_count_by_species(
        sample_count_table_path,
        species_tree.tip_names,
    )
    population_by_branch = _build_species_population_map(
        species_tree,
        node_age_by_id=node_age_by_id,
        default_population_size=default_population_size,
        population_size_table_path=population_size_table_path,
    )
    rng = random.Random(seed)  # nosec B311
    event_rows: list[MultispeciesCoalescentEventRow] = []
    branch_rows: list[MultispeciesCoalescentBranchRow] = []

    def visit(node: TreeNode) -> list[_ActiveGeneLineage]:
        if node.is_leaf():
            lineages = [
                _ActiveGeneLineage(
                    node=TreeNode(name=f"{node.name}__{sample_index}"),
                    age=0.0,
                )
                for sample_index in range(
                    1, sample_count_by_species[node.name or ""] + 1
                )
            ]
        else:
            lineages = []
            for child in node.children:
                lineages.extend(visit(child))
        population = population_by_branch[_species_branch_key(node)]
        return _simulate_gene_lineages_within_population(
            lineages,
            population=population,
            rng=rng,
            event_rows=event_rows,
            branch_rows=branch_rows,
        )

    root_lineages = visit(species_tree.root)
    if len(root_lineages) != 1:
        raise ValueError(
            "multispecies coalescent simulation failed to resolve one gene root"
        )
    gene_root = root_lineages[0].node
    gene_root.branch_length = None
    gene_tree = PhyloTree(root=gene_root, source_format="newick", rooted=True)
    gene_tree.refresh()
    report = MultispeciesCoalescentReport(
        model="multispecies-coalescent",
        species_tree_path=species_tree_path,
        seed=seed,
        species_tip_count=species_tree.tip_count,
        gene_tip_count=gene_tree.tip_count,
        default_population_size=_round_float(default_population_size),
        deep_coalescence_total=sum(
            row.extra_lineage_count
            for row in branch_rows
            if row.included_in_deep_coalescence_total
        ),
        sample_rows=[
            MultispeciesCoalescentSampleRow(
                species_taxon=species_taxon,
                sample_count=sample_count,
                gene_taxa=[
                    f"{species_taxon}__{sample_index}"
                    for sample_index in range(1, sample_count + 1)
                ],
            )
            for species_taxon, sample_count in sorted(sample_count_by_species.items())
        ],
        branch_rows=branch_rows,
        event_rows=event_rows,
    )
    return gene_tree, report


def _simulate_gene_lineages_within_population(
    lineages: list[_ActiveGeneLineage],
    *,
    population: _SpeciesPopulation,
    rng: random.Random,
    event_rows: list[MultispeciesCoalescentEventRow],
    branch_rows: list[MultispeciesCoalescentBranchRow],
) -> list[_ActiveGeneLineage]:
    active_lineages = list(lineages)
    current_age = population.branch_start_age
    entering_count = len(active_lineages)
    event_count_before = len(event_rows)

    while len(active_lineages) > 1:
        coalescent_rate = (
            len(active_lineages)
            * (len(active_lineages) - 1)
            / (2.0 * population.population_size)
        )
        waiting_time = rng.expovariate(coalescent_rate)
        event_age = current_age + waiting_time
        if (
            population.branch_end_age is not None
            and event_age >= population.branch_end_age
        ):
            break
        left_index, right_index = _choose_two_indices(rng, len(active_lineages))
        right = active_lineages.pop(right_index)
        left = active_lineages.pop(left_index)
        left.node.branch_length = _round_float(event_age - left.age)
        right.node.branch_length = _round_float(event_age - right.age)
        parent = TreeNode(children=[left.node, right.node])
        coalesced = _ActiveGeneLineage(node=parent, age=_round_float(event_age))
        active_lineages.append(coalesced)
        current_age = event_age
        event_rows.append(
            MultispeciesCoalescentEventRow(
                event_index=len(event_rows) + 1,
                species_branch=population.species_branch,
                branch_role=population.branch_role,
                descendant_species=population.descendant_species,
                population_size=population.population_size,
                branch_start_age=population.branch_start_age,
                branch_end_age=population.branch_end_age,
                event_age=_round_float(event_age),
                waiting_time=_round_float(waiting_time),
                input_lineage_count=len(active_lineages) + 1,
                output_lineage_count=len(active_lineages),
                left_gene_clade=left.node.descendant_taxa,
                right_gene_clade=right.node.descendant_taxa,
                resulting_gene_clade=coalesced.node.descendant_taxa,
            )
        )

    branch_rows.append(
        MultispeciesCoalescentBranchRow(
            species_branch=population.species_branch,
            branch_role=population.branch_role,
            descendant_species=population.descendant_species,
            branch_duration=population.branch_duration,
            population_size=population.population_size,
            lineage_count_entering=entering_count,
            coalescent_event_count=len(event_rows) - event_count_before,
            lineage_count_exiting=len(active_lineages),
            extra_lineage_count=(
                0
                if not population.included_in_deep_coalescence_total
                else max(len(active_lineages) - 1, 0)
            ),
            included_in_deep_coalescence_total=population.included_in_deep_coalescence_total,
        )
    )
    return active_lineages


def _validate_species_tree(tree: PhyloTree, tree_path: Path) -> None:
    if len(tree.root.children) != 2:
        raise ValueError(
            "multispecies coalescent simulation requires a rooted species tree"
        )
    if any(node.name is None for node in tree.iter_leaves()):
        raise ValueError(
            "multispecies coalescent simulation requires named species-tree tips"
        )
    if len(set(tree.tip_names)) != tree.tip_count:
        raise ValueError(
            "multispecies coalescent simulation requires unique species-tree tip labels"
        )
    if any(
        len(node.children) != 2
        for node in tree.iter_internal_nodes()
        if node is not tree.root
    ):
        raise ValueError(
            "multispecies coalescent simulation requires a strictly binary species tree"
        )
    if any(
        node.branch_length is None
        for node in tree.iter_nodes()
        if node is not tree.root
    ):
        raise ValueError(
            "multispecies coalescent simulation requires complete species-tree branch lengths"
        )
    if not tree_path.exists():
        raise FileNotFoundError(f"tree file not found: {tree_path}")


def _node_age_by_id(tree: PhyloTree, *, root_age: float) -> dict[str, float]:
    node_ages: dict[str, float] = {}

    def visit(node: TreeNode, depth_from_root: float) -> None:
        node_ages[node.node_id or ""] = _round_float(root_age - depth_from_root)
        for child in node.children:
            visit(child, depth_from_root + float(child.branch_length or 0.0))

    visit(tree.root, 0.0)
    return node_ages


def _species_branch_key(node: TreeNode) -> str:
    if node.is_leaf():
        return node.name or "<unnamed>"
    return "|".join(node.descendant_taxa)


def _species_branch_role(node: TreeNode) -> str:
    if node.parent is None:
        return "root-population"
    if node.is_leaf():
        return "tip-branch"
    return "internal-branch"


def _build_species_population_map(
    tree: PhyloTree,
    *,
    node_age_by_id: dict[str, float],
    default_population_size: float,
    population_size_table_path: Path | None,
) -> dict[str, _SpeciesPopulation]:
    override_by_branch = _load_population_size_by_branch(
        population_size_table_path,
        valid_branch_keys=[_species_branch_key(node) for node in tree.iter_nodes()],
    )
    population_map: dict[str, _SpeciesPopulation] = {}
    for node in tree.iter_nodes(order="postorder"):
        species_branch = _species_branch_key(node)
        population_map[species_branch] = _SpeciesPopulation(
            species_branch=species_branch,
            branch_role=_species_branch_role(node),
            descendant_species=node.descendant_taxa,
            population_size=_round_float(
                override_by_branch.get(species_branch, default_population_size)
            ),
            branch_start_age=node_age_by_id[node.node_id or ""],
            branch_end_age=(
                None
                if node.parent is None
                else node_age_by_id[node.parent.node_id or ""]
            ),
        )
    return population_map


def _load_sample_count_by_species(
    sample_count_table_path: Path | None,
    species_tip_names: list[str],
) -> dict[str, int]:
    sample_count_by_species = dict.fromkeys(species_tip_names, 1)
    if sample_count_table_path is None:
        return sample_count_by_species
    table = load_taxon_table(sample_count_table_path)
    if "sample_count" not in table.columns:
        raise ValueError("sample count table must include a 'sample_count' column")
    known_species = set(species_tip_names)
    for row in table.rows:
        species_taxon = row[table.taxon_column]
        if species_taxon not in known_species:
            raise ValueError(
                f"sample count table contains species '{species_taxon}' not present in the species tree"
            )
        try:
            sample_count = int(row["sample_count"])
        except ValueError as error:
            raise ValueError(
                f"sample count for species '{species_taxon}' must be an integer"
            ) from error
        if sample_count < 1:
            raise ValueError(
                f"sample count for species '{species_taxon}' must be at least 1"
            )
        sample_count_by_species[species_taxon] = sample_count
    return sample_count_by_species


def _load_population_size_by_branch(
    population_size_table_path: Path | None,
    *,
    valid_branch_keys: list[str],
) -> dict[str, float]:
    if population_size_table_path is None:
        return {}
    table = load_taxon_table(
        population_size_table_path,
        taxon_column="species_branch",
    )
    if "population_size" not in table.columns:
        raise ValueError(
            "population size table must include a 'population_size' column"
        )
    known_branch_keys = set(valid_branch_keys)
    population_size_by_branch: dict[str, float] = {}
    for row in table.rows:
        species_branch = row["species_branch"]
        if species_branch not in known_branch_keys:
            raise ValueError(
                f"population size table contains unknown species branch '{species_branch}'"
            )
        try:
            population_size = float(row["population_size"])
        except ValueError as error:
            raise ValueError(
                f"population size for species branch '{species_branch}' must be numeric"
            ) from error
        if population_size <= 0.0:
            raise ValueError(
                f"population size for species branch '{species_branch}' must be positive"
            )
        population_size_by_branch[species_branch] = population_size
    return population_size_by_branch


def _choose_two_indices(rng: random.Random, count: int) -> tuple[int, int]:
    left = rng.randrange(count)
    right = rng.randrange(count - 1)
    if right >= left:
        right += 1
    return tuple(sorted((left, right)))
