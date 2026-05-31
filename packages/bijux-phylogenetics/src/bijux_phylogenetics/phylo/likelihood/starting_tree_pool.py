from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodStartingTreePoolReport,
    NucleotideLikelihoodStartingTreeSummary,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.topology import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .multi_start_search import build_random_likelihood_start_tree
from .stepwise_addition import build_likelihood_stepwise_addition_tree
from .starting_tree_validation import (
    validate_nucleotide_likelihood_starting_tree,
)
from .topology_search import (
    evaluate_selected_nucleotide_log_likelihood_from_patterns,
    normalize_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_records,
    resolve_nucleotide_topology_search_surface,
    resolve_nucleotide_topology_search_tree,
)


@dataclass(frozen=True, slots=True)
class _PreparedStartingTree:
    """One generated starting tree with durable source identity."""

    tree_id: str
    source_strategy: str
    generation_seed: int | None
    tree: PhyloTree


def build_nucleotide_likelihood_starting_tree_pool(
    tree: PhyloTree | Path,
    records: list[AlignmentRecord] | Path,
    *,
    model_name: str,
    random_start_tree_count: int = 2,
    random_start_tree_seed: int = 1,
) -> NucleotideLikelihoodStartingTreePoolReport:
    """Build and score one deterministic pool of distinct likelihood starting trees."""
    normalized_model_name = validate_nucleotide_likelihood_starting_tree_pool_model(
        model_name
    )
    validated_random_start_tree_count = (
        validate_nucleotide_likelihood_random_start_tree_count(random_start_tree_count)
    )
    resolved_tree, resolved_tree_path = resolve_nucleotide_topology_search_tree(tree)
    resolved_records, resolved_alignment_path = resolve_nucleotide_topology_search_records(
        records
    )
    normalized_records, compressed_patterns = normalize_nucleotide_topology_search_records(
        resolved_records,
        owner_name="nucleotide likelihood starting-tree pool",
    )
    validate_nucleotide_likelihood_starting_tree(
        resolved_tree,
        compressed_patterns,
        model_name=normalized_model_name,
        workflow_name="nucleotide likelihood starting-tree pool",
    )
    prepared_trees = _prepare_distinct_starting_trees(
        resolved_tree,
        normalized_records,
        compressed_patterns=compressed_patterns,
        model_name=normalized_model_name,
        random_start_tree_count=validated_random_start_tree_count,
        random_start_tree_seed=random_start_tree_seed,
    )
    summaries = [
        _score_prepared_starting_tree(
            prepared_tree,
            normalized_records=normalized_records,
            compressed_patterns=compressed_patterns,
            model_name=normalized_model_name,
        )
        for prepared_tree in prepared_trees
    ]
    return NucleotideLikelihoodStartingTreePoolReport(
        algorithm="nucleotide-likelihood-starting-tree-pool",
        model_name=normalized_model_name.upper(),
        tree_path=None if resolved_tree_path is None else str(resolved_tree_path),
        alignment_path=None
        if resolved_alignment_path is None
        else str(resolved_alignment_path),
        taxon_count=len(compressed_patterns.taxon_order),
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        random_start_tree_count=validated_random_start_tree_count,
        random_start_tree_seed=random_start_tree_seed,
        starting_tree_summaries=summaries,
    )


def build_nucleotide_likelihood_starting_tree_pool_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    random_start_tree_count: int = 2,
    random_start_tree_seed: int = 1,
) -> NucleotideLikelihoodStartingTreePoolReport:
    """Build and score one deterministic likelihood starting-tree pool from paths."""
    return build_nucleotide_likelihood_starting_tree_pool(
        tree_path,
        alignment_path,
        model_name=model_name,
        random_start_tree_count=random_start_tree_count,
        random_start_tree_seed=random_start_tree_seed,
    )


def validate_nucleotide_likelihood_random_start_tree_count(
    random_start_tree_count: int,
) -> int:
    """Require at least one random-seeded start so the pool spans strategies and seeds."""
    if random_start_tree_count < 1:
        raise ValueError("random_start_tree_count must be at least one")
    return random_start_tree_count


def validate_nucleotide_likelihood_starting_tree_pool_model(model_name: str) -> str:
    """Validate the supported likelihood surface for native start-tree pooling."""
    normalized_model_name = model_name.strip().lower()
    if normalized_model_name != "jc69":
        raise ValueError(
            "nucleotide likelihood starting-tree pool model_name must be 'jc69'"
        )
    return normalized_model_name


def _prepare_distinct_starting_trees(
    input_tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    compressed_patterns: CompressedAlignmentSitePatterns,
    model_name: str,
    random_start_tree_count: int,
    random_start_tree_seed: int,
) -> list[_PreparedStartingTree]:
    ordered_taxa = sorted(input_tree.tip_names)
    stepwise_tree, _stepwise_report = build_likelihood_stepwise_addition_tree(
        records,
        model_name="jc69",
    )
    prepared_trees = [
        _PreparedStartingTree(
            tree_id="input-tree",
            source_strategy="input-tree",
            generation_seed=None,
            tree=input_tree.copy().refresh(),
        ),
        _PreparedStartingTree(
            tree_id="likelihood-stepwise-addition-tree",
            source_strategy="likelihood-stepwise-addition-tree",
            generation_seed=None,
            tree=stepwise_tree.refresh(),
        ),
    ]
    for seed in range(
        random_start_tree_seed,
        random_start_tree_seed + random_start_tree_count,
    ):
        prepared_trees.append(
            _PreparedStartingTree(
                tree_id=f"random-tree-seed-{seed}",
                source_strategy="random-tree",
                generation_seed=seed,
                tree=build_random_likelihood_start_tree(ordered_taxa, seed=seed),
            )
        )
    _assert_unique_topology_hashes(
        prepared_trees,
        compressed_patterns=compressed_patterns,
        model_name=model_name,
    )
    return prepared_trees


def _assert_unique_topology_hashes(
    prepared_trees: list[_PreparedStartingTree],
    *,
    compressed_patterns: CompressedAlignmentSitePatterns,
    model_name: str,
) -> None:
    tree_id_by_hash: dict[str, str] = {}
    for prepared_tree in prepared_trees:
        validate_nucleotide_likelihood_starting_tree(
            prepared_tree.tree,
            compressed_patterns,
            model_name=model_name,
            workflow_name="nucleotide likelihood starting-tree pool",
        )
        topology_hash = rooted_topology_fingerprint(prepared_tree.tree)
        other_tree_id = tree_id_by_hash.get(topology_hash)
        if other_tree_id is not None:
            raise ValueError(
                "starting tree pool contains duplicate topology hash "
                f"{topology_hash} for tree_id '{other_tree_id}' and "
                f"tree_id '{prepared_tree.tree_id}'"
            )
        tree_id_by_hash[topology_hash] = prepared_tree.tree_id


def _score_prepared_starting_tree(
    prepared_tree: _PreparedStartingTree,
    *,
    normalized_records: list[AlignmentRecord],
    compressed_patterns: CompressedAlignmentSitePatterns,
    model_name: str,
) -> NucleotideLikelihoodStartingTreeSummary:
    resolved_surface = resolve_nucleotide_topology_search_surface(
        prepared_tree.tree,
        normalized_records,
        model_name=model_name,
    )
    starting_log_likelihood = evaluate_selected_nucleotide_log_likelihood_from_patterns(
        prepared_tree.tree,
        compressed_patterns,
        specification=resolved_surface.specification,
    )
    return NucleotideLikelihoodStartingTreeSummary(
        tree_id=prepared_tree.tree_id,
        source_strategy=prepared_tree.source_strategy,
        generation_seed=prepared_tree.generation_seed,
        topology_hash=rooted_topology_fingerprint(prepared_tree.tree),
        starting_log_likelihood=starting_log_likelihood,
        substitution_parameter_policy=resolved_surface.substitution_parameter_policy,
        substitution_parameter_values=dict(
            resolved_surface.substitution_parameter_values
        ),
        substitution_parameter_warnings=list(
            resolved_surface.substitution_parameter_warnings
        ),
        tree_newick=dumps_newick(prepared_tree.tree),
    )
