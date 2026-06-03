from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .fitch import _leaf_taxa, _resolve_tree
from .matrix import load_parsimony_character_matrix
from .models import (
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
    WagnerCharacterScore,
    WagnerNodeCost,
    WagnerScoreReport,
)
from .weights import resolve_parsimony_character_weights


@dataclass(frozen=True, slots=True)
class _CharacterStateOrder:
    labels: list[str]
    coordinates: dict[str, int]


def score_wagner(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    state_order: list[str] | None = None,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
) -> WagnerScoreReport:
    """Score one discrete character matrix on one tree with ordered Wagner parsimony."""
    resolved_tree, tree_path = _resolve_tree(tree)
    resolved_matrix = (
        matrix
        if isinstance(matrix, FitchCharacterMatrix)
        else load_parsimony_character_matrix(matrix)
    )
    matrix_path = resolved_matrix.matrix_path
    resolved_weights = resolve_parsimony_character_weights(
        resolved_matrix.character_ids,
        character_weights,
    )
    leaf_taxa = _leaf_taxa(resolved_tree)
    missing_from_matrix = sorted(set(leaf_taxa) - set(resolved_matrix.states_by_taxon))
    extra_in_matrix = sorted(set(resolved_matrix.states_by_taxon) - set(leaf_taxa))
    if missing_from_matrix:
        raise ParsimonyAnalysisError(
            "ordered wagner scoring requires every tree taxon to be present in the character matrix",
            code="parsimony_matrix_missing_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "missing_taxa": missing_from_matrix,
            },
        )
    if extra_in_matrix:
        raise ParsimonyAnalysisError(
            "ordered wagner scoring requires matrix taxa to match the tree tips exactly",
            code="parsimony_matrix_extra_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "extra_taxa": extra_in_matrix,
            },
        )

    step_rows: list[WagnerCharacterScore] = []
    node_cost_rows: list[WagnerNodeCost] = []
    total_cost = 0
    total_weighted_score = 0.0
    preorder_nodes = list(resolved_tree.iter_nodes(order="preorder"))
    for character_id in resolved_matrix.character_ids:
        observed_states = sorted(
            {
                states_by_character[character_id]
                for states_by_character in resolved_matrix.states_by_taxon.values()
            }
        )
        resolved_state_order = resolve_wagner_character_state_order(
            character_id=character_id,
            observed_states=observed_states,
            state_order=state_order,
        )
        node_cost_vectors = _score_character(
            resolved_tree.root,
            character_id=character_id,
            states_by_taxon=resolved_matrix.states_by_taxon,
            state_order=resolved_state_order,
        )
        root_key = resolved_tree.root.node_id or node_signature(resolved_tree.root)
        root_cost_vector = node_cost_vectors[root_key]
        character_cost = min(root_cost_vector.values())
        optimal_root_states = [
            state
            for state in resolved_state_order.labels
            if root_cost_vector[state] == character_cost
        ]
        total_cost += character_cost
        character_weight = resolved_weights.weights_by_character[character_id]
        weighted_score = character_cost * character_weight
        total_weighted_score += weighted_score
        step_rows.append(
            WagnerCharacterScore(
                character_id=character_id,
                weighted_step_count=character_cost,
                observed_states=observed_states,
                state_order=resolved_state_order.labels,
                optimal_root_states=optimal_root_states,
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
            for state in resolved_state_order.labels:
                node_cost_rows.append(
                    WagnerNodeCost(
                        character_id=character_id,
                        node=node_signature(node),
                        node_name=node.name,
                        descendant_taxa=node_descendant_taxa(node),
                        state=state,
                        cost=cost_vector[state],
                        is_optimal_state=cost_vector[state] == best_cost,
                    )
                )
    return WagnerScoreReport(
        algorithm="ordered-wagner",
        tree_path=tree_path,
        matrix_path=matrix_path,
        taxon_column=resolved_matrix.taxon_column,
        taxon_count=len(leaf_taxa),
        character_count=resolved_matrix.character_count,
        total_cost=total_cost,
        weights_path=resolved_weights.weights_path,
        total_weighted_score=total_weighted_score,
        step_rows=step_rows,
        node_cost_rows=node_cost_rows,
    )


def resolve_wagner_character_state_order(
    *,
    character_id: str,
    observed_states: list[str],
    state_order: list[str] | None,
) -> _CharacterStateOrder:
    if state_order is not None:
        duplicate_states = sorted(
            {state for state in state_order if state_order.count(state) > 1}
        )
        if duplicate_states:
            raise ParsimonyAnalysisError(
                "ordered wagner scoring requires one unique explicit state order",
                code="parsimony_state_order_invalid",
                details={
                    "character_id": character_id,
                    "duplicate_states": duplicate_states,
                },
            )
        missing_states = sorted(set(observed_states) - set(state_order))
        if missing_states:
            raise ParsimonyAnalysisError(
                "ordered wagner scoring explicit state order must include every observed state",
                code="parsimony_state_order_invalid",
                details={
                    "character_id": character_id,
                    "missing_states": missing_states,
                },
            )
        labels = list(state_order)
        return _CharacterStateOrder(
            labels=labels,
            coordinates={state: index for index, state in enumerate(labels)},
        )
    try:
        ordered_pairs = sorted((int(state), state) for state in observed_states)
    except ValueError as error:
        raise ParsimonyAnalysisError(
            "ordered wagner scoring requires ordinal state labels unless an explicit state order is supplied",
            code="parsimony_state_order_required",
            details={
                "character_id": character_id,
                "observed_states": observed_states,
            },
        ) from error
    labels = [state for _, state in ordered_pairs]
    return _CharacterStateOrder(
        labels=labels,
        coordinates={state: value for value, state in ordered_pairs},
    )


def _score_character(
    root: TreeNode,
    *,
    character_id: str,
    states_by_taxon: dict[str, dict[str, str]],
    state_order: _CharacterStateOrder,
) -> dict[str, dict[str, int]]:
    node_costs: dict[str, dict[str, int]] = {}
    for node in root.iter_nodes(order="postorder"):
        node_key = node.node_id or node_signature(node)
        if node.is_leaf():
            if node.name is None:
                raise ParsimonyAnalysisError(
                    "ordered wagner scoring requires every leaf to have a taxon label",
                    code="parsimony_tree_unnamed_tip",
                )
            observed_state = states_by_taxon[node.name][character_id]
            node_costs[node_key] = {
                state: 0 if state == observed_state else 10**9
                for state in state_order.labels
            }
            continue
        child_cost_vectors = [
            node_costs[child.node_id or node_signature(child)]
            for child in node.children
        ]
        node_costs[node_key] = {}
        for parent_state in state_order.labels:
            parent_coordinate = state_order.coordinates[parent_state]
            total_cost = 0
            for child_cost_vector in child_cost_vectors:
                total_cost += min(
                    child_cost_vector[child_state]
                    + abs(parent_coordinate - state_order.coordinates[child_state])
                    for child_state in state_order.labels
                )
            node_costs[node_key][parent_state] = total_cost
    return node_costs
