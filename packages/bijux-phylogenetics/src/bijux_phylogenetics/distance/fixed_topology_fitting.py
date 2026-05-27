from __future__ import annotations

from dataclasses import dataclass
import math

import numpy

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode, descendant_taxa
from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError

from .fixed_topology_policy import validate_fixed_topology_distance_input

_FIT_ROUND_DIGITS = 12


@dataclass(slots=True)
class FixedTopologyLeastSquaresFit:
    """Least-squares branch-length fit for one fixed tree topology."""

    fitted_tree: PhyloTree
    ordered_taxa: list[str]
    ordered_pairs: list[tuple[str, str]]
    pair_count: int
    branch_count: int
    branch_nodes: list[TreeNode]
    observed_pair_distances: dict[tuple[str, str], float]
    fitted_pair_distances: dict[tuple[str, str], float]
    pair_residuals: dict[tuple[str, str], float]
    pair_weights: dict[tuple[str, str], float]
    residual_sum_squares: float
    objective_sum_squares: float
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


def _validate_pair_weights(
    ordered_pairs: list[tuple[str, str]],
    pair_weights: dict[tuple[str, str], float] | None,
) -> dict[tuple[str, str], float]:
    if pair_weights is None:
        return {pair: 1.0 for pair in ordered_pairs}
    resolved_weights: dict[tuple[str, str], float] = {}
    for pair in ordered_pairs:
        if pair not in pair_weights:
            raise InvalidDistanceMatrixError(
                "least-squares branch fitting requires one weight for every observed taxon pair",
                details={"missing_pair": pair},
            )
        value = float(pair_weights[pair])
        if not math.isfinite(value) or value <= 0.0:
            raise InvalidDistanceMatrixError(
                "least-squares branch fitting requires strictly positive finite pair weights",
                details={"invalid_pair": pair, "invalid_weight": value},
            )
        resolved_weights[pair] = value
    extra_pairs = sorted(set(pair_weights) - set(ordered_pairs))
    if extra_pairs:
        raise InvalidDistanceMatrixError(
            "least-squares branch fitting weights may not include extra taxon pairs",
            details={"extra_pairs": extra_pairs},
        )
    return resolved_weights


def fit_fixed_topology_branch_lengths(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
    *,
    pair_weights: dict[tuple[str, str], float] | None = None,
) -> FixedTopologyLeastSquaresFit:
    """Fit branch lengths for one fixed topology against one full distance matrix."""
    validate_fixed_topology_distance_input(tree, identifiers, distance_lookup)

    branch_nodes = _ordered_branch_nodes(tree)
    ordered_pairs = _pair_order(identifiers)
    resolved_pair_weights = _validate_pair_weights(ordered_pairs, pair_weights)
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
    weight_vector = numpy.array(
        [resolved_pair_weights[pair] for pair in ordered_pairs],
        dtype=float,
    )
    weighted_incidence = incidence * numpy.sqrt(weight_vector)[:, numpy.newaxis]
    weighted_observed = observed * numpy.sqrt(weight_vector)
    solution, _residuals, rank, singular_values = numpy.linalg.lstsq(
        weighted_incidence,
        weighted_observed,
        rcond=None,
    )
    predicted = incidence @ solution
    residuals = observed - predicted
    residual_sum_squares = float(numpy.dot(residuals, residuals))
    objective_sum_squares = float(numpy.dot(weight_vector * residuals, residuals))

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
    observed_pair_distances = {
        pair: round(float(distance_lookup[pair]), _FIT_ROUND_DIGITS)
        for pair in ordered_pairs
    }
    pair_residuals = {
        pair: round(float(residual), _FIT_ROUND_DIGITS)
        for pair, residual in zip(ordered_pairs, residuals, strict=True)
    }
    rounded_pair_weights = {
        pair: round(float(weight), _FIT_ROUND_DIGITS)
        for pair, weight in resolved_pair_weights.items()
    }
    return FixedTopologyLeastSquaresFit(
        fitted_tree=fitted_tree,
        ordered_taxa=list(identifiers),
        ordered_pairs=list(ordered_pairs),
        pair_count=len(ordered_pairs),
        branch_count=len(fitted_branch_nodes),
        branch_nodes=fitted_branch_nodes,
        observed_pair_distances=observed_pair_distances,
        fitted_pair_distances=fitted_pair_distances,
        pair_residuals=pair_residuals,
        pair_weights=rounded_pair_weights,
        residual_sum_squares=round(residual_sum_squares, _FIT_ROUND_DIGITS),
        objective_sum_squares=round(objective_sum_squares, _FIT_ROUND_DIGITS),
        matrix_rank=int(rank),
        singular_values=[
            round(float(value), _FIT_ROUND_DIGITS)
            for value in singular_values.tolist()
        ],
    )
