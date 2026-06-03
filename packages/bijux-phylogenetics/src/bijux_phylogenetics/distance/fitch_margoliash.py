from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, descendant_taxa
from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError

from .fixed_topology_fitting import fit_fixed_topology_branch_lengths
from .imported import (
    _distance_lookup_from_imported,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from .models import FitchMargoliashBranchFit, FitchMargoliashFitReport


def _fitch_margoliash_pair_weights(
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
    *,
    weighting_power: float,
) -> dict[tuple[str, str], float]:
    if not math.isfinite(weighting_power) or weighting_power <= 0.0:
        raise InvalidDistanceMatrixError(
            "fitch-margoliash fitting requires a strictly positive finite weighting power",
            details={"invalid_weighting_power": weighting_power},
        )
    weights: dict[tuple[str, str], float] = {}
    for left_index, left_identifier in enumerate(identifiers):
        for right_identifier in identifiers[left_index + 1 :]:
            observed_distance = float(
                distance_lookup[(left_identifier, right_identifier)]
            )
            if not math.isfinite(observed_distance) or observed_distance <= 0.0:
                raise InvalidDistanceMatrixError(
                    "fitch-margoliash fitting requires strictly positive finite off-diagonal observed distances",
                    details={
                        "invalid_pair": (left_identifier, right_identifier),
                        "invalid_distance": observed_distance,
                    },
                )
            weights[(left_identifier, right_identifier)] = observed_distance ** (
                -weighting_power
            )
    return weights


def fit_fitch_margoliash_tree(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
    *,
    weighting_power: float = 2.0,
) -> tuple[PhyloTree, FitchMargoliashFitReport]:
    """Fit one fixed topology by classical Fitch-Margoliash weighted least squares."""
    fit = fit_fixed_topology_branch_lengths(
        tree,
        identifiers,
        distance_lookup,
        pair_weights=_fitch_margoliash_pair_weights(
            identifiers,
            distance_lookup,
            weighting_power=weighting_power,
        ),
    )
    branch_fits = [
        FitchMargoliashBranchFit(
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
    return fit.fitted_tree, FitchMargoliashFitReport(
        taxa=list(fit.ordered_taxa),
        pair_count=fit.pair_count,
        branch_count=fit.branch_count,
        weighting_power=round(float(weighting_power), 12),
        residual_sum_squares=fit.residual_sum_squares,
        weighted_residual_sum_squares=fit.objective_sum_squares,
        matrix_rank=fit.matrix_rank,
        negative_branch_count=negative_branch_count,
        branch_fits=branch_fits,
    )


def fit_fitch_margoliash_tree_from_imported_distance_matrix(
    matrix_path: Path,
    tree_path: Path,
    *,
    weighting_power: float = 2.0,
) -> tuple[PhyloTree, FitchMargoliashFitReport]:
    """Fit one on-disk tree topology to one imported distance matrix by Fitch-Margoliash."""
    entries = load_imported_distance_matrix(matrix_path)
    validation = validate_imported_distance_matrix(matrix_path)
    distance_lookup, _missing_distance_policy_report = _distance_lookup_from_imported(
        validation,
        entries,
    )
    tree = load_tree(tree_path)
    return fit_fitch_margoliash_tree(
        tree,
        validation.identifiers,
        distance_lookup,
        weighting_power=weighting_power,
    )
