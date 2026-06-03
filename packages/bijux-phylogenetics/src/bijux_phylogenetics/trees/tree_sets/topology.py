from __future__ import annotations

from bijux_phylogenetics.phylo.topology.clades import (
    informative_rooted_clades,
    informative_unrooted_splits,
    robinson_foulds_metrics,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def _format_clade(clade: frozenset[str]) -> str:
    return "|".join(sorted(clade))


def _rooted_topology_id(tree: PhyloTree, shared_taxa: set[str]) -> str:
    return "||".join(
        sorted(
            _format_clade(clade)
            for clade in informative_rooted_clades(tree, shared_taxa)
        )
    )


def _unrooted_topology_id(tree: PhyloTree, shared_taxa: set[str]) -> str:
    return "||".join(
        sorted(
            _format_clade(clade)
            for clade in informative_unrooted_splits(tree, shared_taxa)
        )
    )


def _clade_signature(tree: PhyloTree, shared_taxa: set[str], taxon: str) -> str:
    containing_clades = sorted(
        _format_clade(clade)
        for clade in informative_rooted_clades(tree, shared_taxa)
        if taxon in clade
    )
    if not containing_clades:
        return "(singleton-placement)"
    return "||".join(containing_clades)


def _clade_counts(
    trees: list[PhyloTree], shared_taxa: set[str]
) -> dict[frozenset[str], int]:
    counts: dict[frozenset[str], int] = {}
    for tree in trees:
        for clade in informative_rooted_clades(tree, shared_taxa):
            counts[clade] = counts.get(clade, 0) + 1
    return counts


def _clades_conflict(left: frozenset[str], right: frozenset[str]) -> bool:
    return bool(left & right) and not (left <= right or right <= left)


def _tree_distance(
    left: PhyloTree, right: PhyloTree, shared_taxa: set[str]
) -> tuple[int, float]:
    metrics = robinson_foulds_metrics(
        left,
        right,
        shared_taxa,
        rf_mode="rooted",
    )
    return metrics.distance, metrics.normalized_distance
