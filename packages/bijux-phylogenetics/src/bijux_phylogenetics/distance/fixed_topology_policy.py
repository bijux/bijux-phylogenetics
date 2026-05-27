from __future__ import annotations

from collections import Counter

from bijux_phylogenetics.phylo.topology.distance_joining import (
    validate_distance_lookup,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    DuplicateTaxonError,
    InvalidDistanceMatrixError,
    UnnamedTipError,
)


def validate_fixed_topology_distance_input(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> None:
    """Require named unique tree tips and one exact taxon match to the distance matrix."""
    validate_distance_lookup(identifiers, distance_lookup)
    tip_names = [node.name for node in tree.iter_leaves()]
    unnamed_tip_count = sum(1 for name in tip_names if not name)
    if unnamed_tip_count:
        raise UnnamedTipError(
            f"fixed-topology distance analysis requires named tips and found {unnamed_tip_count} unnamed tips"
        )
    ordered_tip_names = [name for name in tip_names if name is not None]
    duplicate_tip_names = sorted(
        name for name, count in Counter(ordered_tip_names).items() if count > 1
    )
    if duplicate_tip_names:
        raise DuplicateTaxonError(
            "fixed-topology distance analysis requires unique tip labels and found duplicates: "
            + ", ".join(duplicate_tip_names)
        )
    tree_taxon_set = set(tree.tip_names)
    matrix_taxon_set = set(identifiers)
    if tree_taxon_set != matrix_taxon_set:
        raise InvalidDistanceMatrixError(
            "distance matrix taxa do not match the fixed tree tip labels",
            details={
                "tree_only_taxa": sorted(tree_taxon_set - matrix_taxon_set),
                "matrix_only_taxa": sorted(matrix_taxon_set - tree_taxon_set),
            },
        )
