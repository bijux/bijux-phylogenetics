from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, descendant_taxa

from .fixed_topology_fitting import fit_fixed_topology_branch_lengths
from .imported import (
    _distance_lookup_from_imported,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from .models import OrdinaryLeastSquaresBranchFit, OrdinaryLeastSquaresFitReport


def _build_symmetric_matrix(
    identifiers: list[str],
    pair_values: dict[tuple[str, str], float],
) -> list[list[float]]:
    matrix: list[list[float]] = []
    for left_index, left_identifier in enumerate(identifiers):
        row: list[float] = []
        for right_index, right_identifier in enumerate(identifiers):
            if left_index == right_index:
                row.append(0.0)
            elif left_index < right_index:
                row.append(pair_values[(left_identifier, right_identifier)])
            else:
                row.append(pair_values[(right_identifier, left_identifier)])
        matrix.append(row)
    return matrix


def _condition_number(singular_values: list[float]) -> float:
    if not singular_values:
        return math.inf
    smallest = min(singular_values)
    if smallest == 0.0:
        return math.inf
    return round(singular_values[0] / smallest, 12)


def fit_ordinary_least_squares_tree(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, OrdinaryLeastSquaresFitReport]:
    """Fit one fixed topology by ordinary least squares over the path-incidence matrix."""
    fit = fit_fixed_topology_branch_lengths(tree, identifiers, distance_lookup)
    branch_fits = [
        OrdinaryLeastSquaresBranchFit(
            branch_id=node.node_id or "",
            child_name=node.name,
            descendant_taxa=descendant_taxa(node),
            fitted_branch_length=node.branch_length or 0.0,
        )
        for node in fit.branch_nodes
    ]
    negative_branch_count = sum(
        1 for row in branch_fits if row.fitted_branch_length < 0.0
    )
    return fit.fitted_tree, OrdinaryLeastSquaresFitReport(
        taxa=list(fit.ordered_taxa),
        pair_count=fit.pair_count,
        branch_count=fit.branch_count,
        residual_sum_squares=fit.residual_sum_squares,
        matrix_rank=fit.matrix_rank,
        condition_number=_condition_number(fit.singular_values),
        negative_branch_count=negative_branch_count,
        fitted_distance_matrix=_build_symmetric_matrix(
            fit.ordered_taxa,
            fit.fitted_pair_distances,
        ),
        residual_matrix=_build_symmetric_matrix(
            fit.ordered_taxa,
            fit.pair_residuals,
        ),
        branch_fits=branch_fits,
    )


def fit_ordinary_least_squares_tree_from_imported_distance_matrix(
    matrix_path: Path,
    tree_path: Path,
) -> tuple[PhyloTree, OrdinaryLeastSquaresFitReport]:
    """Fit one on-disk tree topology to one imported distance matrix by ordinary least squares."""
    entries = load_imported_distance_matrix(matrix_path)
    validation = validate_imported_distance_matrix(matrix_path)
    distance_lookup, _missing_distance_policy_report = _distance_lookup_from_imported(
        validation,
        entries,
    )
    tree = load_tree(tree_path)
    return fit_ordinary_least_squares_tree(
        tree,
        validation.identifiers,
        distance_lookup,
    )
