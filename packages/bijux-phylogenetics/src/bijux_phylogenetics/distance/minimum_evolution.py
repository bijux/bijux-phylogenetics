from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, descendant_taxa

from .fixed_topology_fitting import fit_fixed_topology_branch_lengths
from .imported import (
    _distance_lookup_from_imported,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from .models import MinimumEvolutionBranchFit, MinimumEvolutionScoreReport


def fit_minimum_evolution_tree(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> tuple[PhyloTree, MinimumEvolutionScoreReport]:
    """Fit one fixed topology to a distance matrix and score it by total branch length."""
    fit = fit_fixed_topology_branch_lengths(tree, identifiers, distance_lookup)
    branch_fits = [
        MinimumEvolutionBranchFit(
            branch_id=node.node_id or "",
            child_name=node.name,
            descendant_taxa=descendant_taxa(node),
            fitted_branch_length=node.branch_length or 0.0,
        )
        for node in fit.branch_nodes
    ]
    total_fitted_branch_length = round(fit.fitted_tree.total_branch_length(), 15)
    negative_branch_count = sum(
        1 for row in branch_fits if row.fitted_branch_length < 0.0
    )
    return fit.fitted_tree, MinimumEvolutionScoreReport(
        taxa=list(fit.ordered_taxa),
        pair_count=fit.pair_count,
        branch_count=fit.branch_count,
        minimum_evolution_score=total_fitted_branch_length,
        total_fitted_branch_length=total_fitted_branch_length,
        negative_branch_count=negative_branch_count,
        branch_fits=branch_fits,
    )


def fit_minimum_evolution_tree_from_imported_distance_matrix(
    matrix_path: Path,
    tree_path: Path,
) -> tuple[PhyloTree, MinimumEvolutionScoreReport]:
    """Fit one on-disk tree topology to one imported distance matrix by minimum evolution."""
    entries = load_imported_distance_matrix(matrix_path)
    validation = validate_imported_distance_matrix(matrix_path)
    distance_lookup, _missing_distance_policy_report = _distance_lookup_from_imported(
        validation,
        entries,
    )
    tree = load_tree(tree_path)
    return fit_minimum_evolution_tree(
        tree,
        validation.identifiers,
        distance_lookup,
    )
