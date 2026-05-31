from __future__ import annotations

from pathlib import Path

import numpy

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, descendant_taxa

from .fixed_topology_fitting import (
    LeastSquaresSystem,
    build_fixed_topology_least_squares_system,
)
from .imported import (
    _distance_lookup_from_imported,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from .models import (
    NonnegativeLeastSquaresActiveConstraint,
    NonnegativeLeastSquaresBranchFit,
    NonnegativeLeastSquaresFitReport,
)
from .ordinary_least_squares import _build_symmetric_matrix, _condition_number

_NNLS_TOLERANCE = 1e-10


def _solve_active_set(
    system: LeastSquaresSystem,
) -> tuple[numpy.ndarray, numpy.ndarray, list[int]]:
    weighted_incidence = (
        system.incidence * numpy.sqrt(system.weight_vector)[:, numpy.newaxis]
    )
    weighted_observed = system.observed * numpy.sqrt(system.weight_vector)
    branch_count = weighted_incidence.shape[1]
    passive: set[int] = set()
    active: set[int] = set(range(branch_count))
    solution = numpy.zeros(branch_count, dtype=float)
    dual = weighted_incidence.T @ (weighted_observed - (weighted_incidence @ solution))

    while (
        active
        and float(numpy.max([dual[index] for index in active], initial=0.0))
        > _NNLS_TOLERANCE
    ):
        entering = max(active, key=lambda index: float(dual[index]))
        passive.add(entering)
        active.remove(entering)

        while True:
            passive_order = sorted(passive)
            passive_matrix = weighted_incidence[:, passive_order]
            passive_solution, _residuals, _rank, _singular_values = numpy.linalg.lstsq(
                passive_matrix,
                weighted_observed,
                rcond=None,
            )
            candidate = numpy.zeros(branch_count, dtype=float)
            for local_index, branch_index in enumerate(passive_order):
                candidate[branch_index] = passive_solution[local_index]
            if all(candidate[index] > _NNLS_TOLERANCE for index in passive):
                solution = candidate
                break

            leaving = [
                index for index in passive if candidate[index] <= _NNLS_TOLERANCE
            ]
            alpha = min(
                solution[index] / (solution[index] - candidate[index])
                for index in leaving
                if solution[index] - candidate[index] > 0.0
            )
            solution = solution + alpha * (candidate - solution)
            exiting = [index for index in passive if solution[index] <= _NNLS_TOLERANCE]
            for index in exiting:
                solution[index] = 0.0
                passive.remove(index)
                active.add(index)
        dual = weighted_incidence.T @ (
            weighted_observed - (weighted_incidence @ solution)
        )
    return solution, system.incidence @ solution, sorted(active)


def fit_nonnegative_least_squares_tree(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, NonnegativeLeastSquaresFitReport]:
    """Fit one fixed topology with nonnegative branch constraints."""
    branch_nodes, system = build_fixed_topology_least_squares_system(
        tree,
        identifiers,
        distance_lookup,
    )
    solution, predicted, active_constraints = _solve_active_set(system)
    residuals = system.observed - predicted
    residual_sum_squares = round(float(numpy.dot(residuals, residuals)), 12)

    fitted_tree = tree.copy()
    fitted_branch_nodes = [
        node
        for node in fitted_tree.iter_nodes(order="preorder")
        if node is not fitted_tree.root
    ]
    for node, branch_length in zip(fitted_branch_nodes, solution, strict=True):
        rounded_branch_length = round(float(branch_length), 12)
        if abs(rounded_branch_length) <= 1e-12:
            rounded_branch_length = 0.0
        node.branch_length = rounded_branch_length
    fitted_tree.refresh()

    branch_fits = [
        NonnegativeLeastSquaresBranchFit(
            branch_id=node.node_id or "",
            child_name=node.name,
            descendant_taxa=descendant_taxa(node),
            fitted_branch_length=node.branch_length or 0.0,
        )
        for node in fitted_branch_nodes
    ]
    active_constraint_rows = [
        NonnegativeLeastSquaresActiveConstraint(
            branch_id=fitted_branch_nodes[index].node_id or "",
            child_name=fitted_branch_nodes[index].name,
            descendant_taxa=descendant_taxa(fitted_branch_nodes[index]),
        )
        for index in active_constraints
    ]
    singular_values = numpy.linalg.svd(system.incidence, compute_uv=False).tolist()
    return fitted_tree, NonnegativeLeastSquaresFitReport(
        taxa=list(identifiers),
        pair_count=len(system.ordered_pairs),
        branch_count=len(fitted_branch_nodes),
        residual_sum_squares=residual_sum_squares,
        condition_number=_condition_number(
            [round(float(value), 12) for value in singular_values]
        ),
        active_constraint_count=len(active_constraint_rows),
        fitted_distance_matrix=_build_symmetric_matrix(
            identifiers,
            {
                pair: round(float(distance), 12)
                for pair, distance in zip(system.ordered_pairs, predicted, strict=True)
            },
        ),
        residual_matrix=_build_symmetric_matrix(
            identifiers,
            {
                pair: round(float(residual), 12)
                for pair, residual in zip(system.ordered_pairs, residuals, strict=True)
            },
        ),
        branch_fits=branch_fits,
        active_constraints=active_constraint_rows,
    )


def fit_nonnegative_least_squares_tree_from_imported_distance_matrix(
    matrix_path: Path,
    tree_path: Path,
) -> tuple[PhyloTree, NonnegativeLeastSquaresFitReport]:
    """Fit one on-disk tree topology to one imported distance matrix with NNLS constraints."""
    entries = load_imported_distance_matrix(matrix_path)
    validation = validate_imported_distance_matrix(matrix_path)
    distance_lookup, _missing_distance_policy_report = _distance_lookup_from_imported(
        validation,
        entries,
    )
    tree = load_tree(tree_path)
    return fit_nonnegative_least_squares_tree(
        tree,
        validation.identifiers,
        distance_lookup,
    )
