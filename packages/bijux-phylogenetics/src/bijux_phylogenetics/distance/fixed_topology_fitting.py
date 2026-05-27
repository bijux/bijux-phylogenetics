from __future__ import annotations

from dataclasses import dataclass

import numpy

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode, descendant_taxa

from .fixed_topology_policy import validate_fixed_topology_distance_input

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
    validate_fixed_topology_distance_input(tree, identifiers, distance_lookup)

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
