from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import numpy

from bijux_phylogenetics.phylo.topology.distance_joining import (
    validate_distance_lookup,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode, descendant_taxa
from bijux_phylogenetics.runtime.errors import (
    DuplicateTaxonError,
    InvalidDistanceMatrixError,
    UnnamedTipError,
)

_FIT_ROUND_DIGITS = 12


@dataclass(slots=True)
class FixedTopologyLeastSquaresFit:
    """Least-squares branch-length fit for one fixed tree topology."""

    fitted_tree: PhyloTree
    ordered_taxa: list[str]
    pair_count: int
    branch_count: int
    branch_nodes: list[TreeNode]
    fitted_pair_distances: dict[tuple[str, str], float]
    residual_sum_squares: float
    matrix_rank: int
    singular_values: list[float]


def _ordered_branch_nodes(tree: PhyloTree) -> list[TreeNode]:
    return [node for node in tree.iter_nodes(order="preorder") if node is not tree.root]


def _validate_tree_tip_labels(tree: PhyloTree) -> None:
    tip_names = [node.name for node in tree.iter_leaves()]
    unnamed_tip_count = sum(1 for name in tip_names if not name)
    if unnamed_tip_count:
        raise UnnamedTipError(
            f"fixed-topology fitting requires named tips and found {unnamed_tip_count} unnamed tips"
        )
    ordered_tip_names = [name for name in tip_names if name is not None]
    duplicate_tip_names = sorted(
        name for name, count in Counter(ordered_tip_names).items() if count > 1
    )
    if duplicate_tip_names:
        raise DuplicateTaxonError(
            "fixed-topology fitting requires unique tip labels and found duplicates: "
            + ", ".join(duplicate_tip_names)
        )


def _validate_tree_taxa(tree: PhyloTree, identifiers: list[str]) -> None:
    _validate_tree_tip_labels(tree)
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


def _pair_order(identifiers: list[str]) -> list[tuple[str, str]]:
    ordered_pairs: list[tuple[str, str]] = []
    for left_index, left_identifier in enumerate(identifiers):
        for right_identifier in identifiers[left_index + 1 :]:
            ordered_pairs.append((left_identifier, right_identifier))
    return ordered_pairs


def fit_fixed_topology_branch_lengths(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> FixedTopologyLeastSquaresFit:
    """Fit branch lengths for one fixed topology against one full distance matrix."""
    validate_distance_lookup(identifiers, distance_lookup)
    _validate_tree_taxa(tree, identifiers)

    branch_nodes = _ordered_branch_nodes(tree)
    ordered_pairs = _pair_order(identifiers)
    descendant_sets = [set(descendant_taxa(node)) for node in branch_nodes]

    incidence_rows: list[list[float]] = []
    observed_distances: list[float] = []
    for left_identifier, right_identifier in ordered_pairs:
        incidence_rows.append(
            [
                1.0
                if ((left_identifier in descendant_set) != (right_identifier in descendant_set))
                else 0.0
                for descendant_set in descendant_sets
            ]
        )
        observed_distances.append(distance_lookup[(left_identifier, right_identifier)])

    incidence = numpy.array(incidence_rows, dtype=float)
    observed = numpy.array(observed_distances, dtype=float)
    solution, _residuals, rank, singular_values = numpy.linalg.lstsq(
        incidence,
        observed,
        rcond=None,
    )
    predicted = incidence @ solution
    residual_sum_squares = float(numpy.dot(observed - predicted, observed - predicted))

    fitted_tree = tree.copy()
    fitted_branch_nodes = _ordered_branch_nodes(fitted_tree)
    for node, branch_length in zip(fitted_branch_nodes, solution, strict=True):
        rounded_branch_length = round(float(branch_length), _FIT_ROUND_DIGITS)
        if abs(rounded_branch_length) <= 1e-12:
            rounded_branch_length = 0.0
        node.branch_length = rounded_branch_length
    fitted_tree.refresh()

    fitted_pair_distances = {
        pair: round(float(distance), _FIT_ROUND_DIGITS)
        for pair, distance in zip(ordered_pairs, predicted, strict=True)
    }
    return FixedTopologyLeastSquaresFit(
        fitted_tree=fitted_tree,
        ordered_taxa=list(identifiers),
        pair_count=len(ordered_pairs),
        branch_count=len(fitted_branch_nodes),
        branch_nodes=fitted_branch_nodes,
        fitted_pair_distances=fitted_pair_distances,
        residual_sum_squares=round(residual_sum_squares, _FIT_ROUND_DIGITS),
        matrix_rank=int(rank),
        singular_values=[
            round(float(value), _FIT_ROUND_DIGITS)
            for value in singular_values.tolist()
        ],
    )
