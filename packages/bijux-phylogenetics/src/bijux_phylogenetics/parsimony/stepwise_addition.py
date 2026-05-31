from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from bijux_phylogenetics.phylo.topology import (
    PhyloTree,
    StepwiseAdditionTreeReport,
    build_greedy_stepwise_addition_tree,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .models import FitchCharacterMatrix, ParsimonyCharacterWeights, SankoffCostMatrix
from .topology_search import (
    resolve_topology_search_cost_matrix,
    resolve_topology_search_matrix,
    resolve_topology_search_method,
    resolve_topology_search_weights,
    score_topology_search_tree,
)


def build_parsimony_stepwise_addition_tree(
    matrix: FitchCharacterMatrix | Path,
    *,
    method: str = "fitch",
    taxon_column: str | None = None,
    insertion_order: list[str] | None = None,
    state_order: list[str] | None = None,
    cost_matrix: SankoffCostMatrix | Path | None = None,
    allow_asymmetric_costs: bool = False,
    character_weights: ParsimonyCharacterWeights
    | Mapping[str, float]
    | Path
    | None = None,
) -> tuple[PhyloTree, StepwiseAdditionTreeReport]:
    """Build one rooted tree by greedy stepwise addition under governed parsimony scoring."""
    workflow_name = "parsimony stepwise addition"
    resolved_matrix = resolve_topology_search_matrix(
        matrix,
        taxon_column=taxon_column,
    )
    resolved_method = resolve_topology_search_method(
        method,
        workflow_name=workflow_name,
    )
    resolved_cost_matrix = resolve_topology_search_cost_matrix(
        method=resolved_method,
        cost_matrix=cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        workflow_name=workflow_name,
    )
    resolved_character_weights = resolve_topology_search_weights(
        resolved_matrix,
        character_weights,
    )
    resolved_insertion_order = _resolve_stepwise_insertion_order(
        resolved_matrix,
        insertion_order=insertion_order,
    )

    def score_tree(tree: PhyloTree) -> float:
        restricted_matrix = _restrict_matrix_to_tree_taxa(
            resolved_matrix, tree.tip_names
        )
        return score_topology_search_tree(
            tree,
            restricted_matrix,
            method=resolved_method,
            state_order=state_order,
            cost_matrix=resolved_cost_matrix,
            allow_asymmetric_costs=allow_asymmetric_costs,
            character_weights=resolved_character_weights,
        )

    return build_greedy_stepwise_addition_tree(
        resolved_insertion_order,
        score_tree=score_tree,
        objective_name=f"parsimony-{resolved_method}",
        objective_direction="minimize",
    )


def _restrict_matrix_to_tree_taxa(
    matrix: FitchCharacterMatrix,
    tip_names: list[str],
) -> FitchCharacterMatrix:
    return FitchCharacterMatrix(
        matrix_path=matrix.matrix_path,
        taxon_column=matrix.taxon_column,
        character_ids=matrix.character_ids,
        states_by_taxon={taxon: matrix.states_by_taxon[taxon] for taxon in tip_names},
    )


def _resolve_stepwise_insertion_order(
    matrix: FitchCharacterMatrix,
    *,
    insertion_order: list[str] | None,
) -> list[str]:
    default_insertion_order = list(matrix.states_by_taxon)
    if insertion_order is None:
        return default_insertion_order

    duplicate_taxa = sorted(
        {taxon for taxon in insertion_order if insertion_order.count(taxon) > 1}
    )
    if duplicate_taxa:
        raise ParsimonyAnalysisError(
            "parsimony stepwise addition requires distinct taxa in insertion order",
            code="parsimony_stepwise_addition_insertion_order_duplicate_taxa",
            details={"duplicate_taxa": duplicate_taxa},
        )
    matrix_taxa = set(default_insertion_order)
    insertion_taxa = set(insertion_order)
    missing_taxa = sorted(matrix_taxa - insertion_taxa)
    unexpected_taxa = sorted(insertion_taxa - matrix_taxa)
    if missing_taxa or unexpected_taxa:
        raise ParsimonyAnalysisError(
            "parsimony stepwise addition insertion order must match matrix taxa exactly",
            code="parsimony_stepwise_addition_insertion_order_taxa_mismatch",
            details={
                "missing_taxa": missing_taxa,
                "unexpected_taxa": unexpected_taxa,
            },
        )
    return list(insertion_order)
