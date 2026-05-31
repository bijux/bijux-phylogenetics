from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .cost_matrix import load_sankoff_cost_matrix, validate_sankoff_cost_matrix
from .fitch import _leaf_taxa, _resolve_tree
from .matrix import load_parsimony_character_matrix
from .models import (
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
    SankoffCharacterScore,
    SankoffCostMatrix,
    SankoffNodeCost,
    SankoffNodeSelection,
    SankoffScoreReport,
)
from .weights import resolve_parsimony_character_weights


def score_sankoff(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    cost_matrix: SankoffCostMatrix | Path,
    *,
    allow_asymmetric_costs: bool = False,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
) -> SankoffScoreReport:
    """Score one discrete character matrix on one tree with Sankoff parsimony."""
    resolved_tree, tree_path = _resolve_tree(tree)
    resolved_matrix = (
        matrix
        if isinstance(matrix, FitchCharacterMatrix)
        else load_parsimony_character_matrix(matrix)
    )
    resolved_cost_matrix = (
        cost_matrix
        if isinstance(cost_matrix, SankoffCostMatrix)
        else load_sankoff_cost_matrix(
            cost_matrix,
            observed_states=_observed_states_across_matrix(resolved_matrix),
            allow_asymmetric_costs=allow_asymmetric_costs,
        )
    )
    if isinstance(cost_matrix, SankoffCostMatrix):
        resolved_cost_matrix = validate_sankoff_cost_matrix(
            resolved_cost_matrix,
            observed_states=_observed_states_across_matrix(resolved_matrix),
            allow_asymmetric_costs=allow_asymmetric_costs,
        )
    matrix_path = resolved_matrix.matrix_path
    cost_matrix_path = resolved_cost_matrix.matrix_path
    resolved_weights = resolve_parsimony_character_weights(
        resolved_matrix.character_ids,
        character_weights,
    )
    leaf_taxa = _leaf_taxa(resolved_tree)
    missing_from_matrix = sorted(set(leaf_taxa) - set(resolved_matrix.states_by_taxon))
    extra_in_matrix = sorted(set(resolved_matrix.states_by_taxon) - set(leaf_taxa))
    if missing_from_matrix:
        raise ParsimonyAnalysisError(
            "sankoff scoring requires every tree taxon to be present in the character matrix",
            code="parsimony_matrix_missing_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "missing_taxa": missing_from_matrix,
            },
        )
    if extra_in_matrix:
        raise ParsimonyAnalysisError(
            "sankoff scoring requires matrix taxa to match the tree tips exactly",
            code="parsimony_matrix_extra_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "extra_taxa": extra_in_matrix,
            },
        )

    step_rows: list[SankoffCharacterScore] = []
    node_cost_rows: list[SankoffNodeCost] = []
    node_selection_rows: list[SankoffNodeSelection] = []
    total_cost = 0.0
    total_weighted_score = 0.0
    preorder_nodes = list(resolved_tree.iter_nodes(order="preorder"))
    for character_id in resolved_matrix.character_ids:
        observed_states = sorted(
            {
                states_by_character[character_id]
                for states_by_character in resolved_matrix.states_by_taxon.values()
            }
        )
        node_cost_vectors = _score_character(
            resolved_tree.root,
            character_id=character_id,
            states_by_taxon=resolved_matrix.states_by_taxon,
            cost_matrix=resolved_cost_matrix,
        )
        root_key = resolved_tree.root.node_id or node_signature(resolved_tree.root)
        root_cost_vector = node_cost_vectors[root_key]
        character_cost = min(root_cost_vector.values())
        total_cost += character_cost
        character_weight = resolved_weights.weights_by_character[character_id]
        weighted_score = character_cost * character_weight
        total_weighted_score += weighted_score
        step_rows.append(
            SankoffCharacterScore(
                character_id=character_id,
                minimum_cost=character_cost,
                observed_states=observed_states,
                matrix_states=resolved_cost_matrix.states,
                character_weight=character_weight,
                weighted_score=weighted_score,
            )
        )
        for node in preorder_nodes:
            if node.is_leaf():
                continue
            node_key = node.node_id or node_signature(node)
            cost_vector = node_cost_vectors[node_key]
            best_cost = min(cost_vector.values())
            optimal_states = [
                state
                for state in resolved_cost_matrix.states
                if cost_vector[state] == best_cost
            ]
            node_selection_rows.append(
                SankoffNodeSelection(
                    character_id=character_id,
                    node=node_signature(node),
                    node_name=node.name,
                    descendant_taxa=node_descendant_taxa(node),
                    optimal_states=optimal_states,
                    tie_states=optimal_states if len(optimal_states) > 1 else [],
                )
            )
            for state in resolved_cost_matrix.states:
                node_cost_rows.append(
                    SankoffNodeCost(
                        character_id=character_id,
                        node=node_signature(node),
                        node_name=node.name,
                        descendant_taxa=node_descendant_taxa(node),
                        state=state,
                        cost=cost_vector[state],
                        is_optimal_state=cost_vector[state] == best_cost,
                    )
                )
    return SankoffScoreReport(
        algorithm="sankoff",
        tree_path=tree_path,
        matrix_path=matrix_path,
        cost_matrix_path=cost_matrix_path,
        taxon_column=resolved_matrix.taxon_column,
        taxon_count=len(leaf_taxa),
        character_count=resolved_matrix.character_count,
        total_cost=total_cost,
        weights_path=resolved_weights.weights_path,
        total_weighted_score=total_weighted_score,
        validation_warnings=resolved_cost_matrix.validation_warnings,
        step_rows=step_rows,
        node_cost_rows=node_cost_rows,
        node_selection_rows=node_selection_rows,
    )


def _observed_states_across_matrix(matrix: FitchCharacterMatrix) -> list[str]:
    return sorted(
        {
            state
            for states_by_character in matrix.states_by_taxon.values()
            for state in states_by_character.values()
        }
    )


def _score_character(
    root: TreeNode,
    *,
    character_id: str,
    states_by_taxon: dict[str, dict[str, str]],
    cost_matrix: SankoffCostMatrix,
) -> dict[str, dict[str, float]]:
    node_costs: dict[str, dict[str, float]] = {}
    impossible_cost = float("inf")
    for node in root.iter_nodes(order="postorder"):
        node_key = node.node_id or node_signature(node)
        if node.is_leaf():
            if node.name is None:
                raise ParsimonyAnalysisError(
                    "sankoff scoring requires every leaf to have a taxon label",
                    code="parsimony_tree_unnamed_tip",
                )
            observed_state = states_by_taxon[node.name][character_id]
            node_costs[node_key] = {
                state: 0.0 if state == observed_state else impossible_cost
                for state in cost_matrix.states
            }
            continue
        child_cost_vectors = [
            node_costs[child.node_id or node_signature(child)]
            for child in node.children
        ]
        node_costs[node_key] = {}
        for parent_state in cost_matrix.states:
            total_cost = 0.0
            for child_cost_vector in child_cost_vectors:
                total_cost += min(
                    child_cost_vector[child_state]
                    + cost_matrix.costs[parent_state][child_state]
                    for child_state in cost_matrix.states
                )
            node_costs[node_key][parent_state] = total_cost
    return node_costs
