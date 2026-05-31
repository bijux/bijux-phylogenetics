from __future__ import annotations

from bijux_phylogenetics.phylo.topology.random_bifurcating import (
    generate_random_bifurcating_tree,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def build_random_likelihood_start_tree(
    ordered_taxa: list[str],
    *,
    seed: int,
) -> PhyloTree:
    """Generate one rooted random start tree and relabel its tips to the target taxa."""
    random_tree, _report = generate_random_bifurcating_tree(
        ordered_taxa,
        seed=seed,
        branch_length_policy="uniform",
    )
    return random_tree
