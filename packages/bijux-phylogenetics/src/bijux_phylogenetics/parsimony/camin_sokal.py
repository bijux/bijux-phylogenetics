from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .fitch import _leaf_taxa, _resolve_tree
from .matrix import load_parsimony_character_matrix
from .models import (
    CaminSokalBranchChange,
    CaminSokalCharacterScore,
    CaminSokalScoreReport,
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
)
from .weights import resolve_parsimony_character_weights

_BINARY_STATES = frozenset({"0", "1"})
_ROOT_STATE = "0"
_STATES = (_ROOT_STATE, "1")
_INFINITE_COST = float("inf")


def score_camin_sokal(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
) -> CaminSokalScoreReport:
    """Score one binary character matrix on one rooted tree with irreversible Camin-Sokal parsimony."""
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
            "camin-sokal scoring requires every tree taxon to be present in the character matrix",
            code="parsimony_matrix_missing_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "missing_taxa": missing_from_matrix,
            },
        )
    if extra_in_matrix:
        raise ParsimonyAnalysisError(
            "camin-sokal scoring requires matrix taxa to match the tree tips exactly",
            code="parsimony_matrix_extra_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "extra_taxa": extra_in_matrix,
            },
        )

    step_rows: list[CaminSokalCharacterScore] = []
    branch_change_rows: list[CaminSokalBranchChange] = []
    total_gains = 0
    total_weighted_score = 0.0
    for character_id in resolved_matrix.character_ids:
        character_weight = resolved_weights.weights_by_character[character_id]
        observed_states = {
            states_by_character[character_id]
            for states_by_character in resolved_matrix.states_by_taxon.values()
        }
        invalid_states = sorted(observed_states - _BINARY_STATES)
        if invalid_states:
            raise ParsimonyAnalysisError(
                "camin-sokal scoring requires binary 0/1 character states unless the trait is binarized first",
                code="parsimony_matrix_multistate_not_binarized",
                details={
                    "character_id": character_id,
                    "invalid_states": invalid_states,
                },
            )
        branch_changes = _score_irreversible_character(
            resolved_tree.root,
            character_id=character_id,
            states_by_taxon=resolved_matrix.states_by_taxon,
        )
        derived_taxon_count = sum(
            1
            for states_by_character in resolved_matrix.states_by_taxon.values()
            if states_by_character[character_id] == "1"
        )
        total_gains += len(branch_changes)
        weighted_score = len(branch_changes) * character_weight
        total_weighted_score += weighted_score
        branch_change_rows.extend(branch_changes)
        step_rows.append(
            CaminSokalCharacterScore(
                character_id=character_id,
                derived_taxon_count=derived_taxon_count,
                gain_count=len(branch_changes),
                root_state=_ROOT_STATE,
                character_weight=character_weight,
                weighted_score=weighted_score,
            )
        )
    return CaminSokalScoreReport(
        algorithm="camin-sokal",
        tree_path=tree_path,
        matrix_path=matrix_path,
        taxon_column=resolved_matrix.taxon_column,
        taxon_count=len(leaf_taxa),
        character_count=resolved_matrix.character_count,
        root_state=_ROOT_STATE,
        total_gains=total_gains,
        weights_path=resolved_weights.weights_path,
        total_weighted_score=total_weighted_score,
        step_rows=step_rows,
        branch_change_rows=branch_change_rows,
    )


def _score_irreversible_character(
    root: TreeNode,
    *,
    character_id: str,
    states_by_taxon: dict[str, dict[str, str]],
) -> list[CaminSokalBranchChange]:
    cost_rows = _compute_irreversible_costs(
        root,
        character_id=character_id,
        states_by_taxon=states_by_taxon,
    )
    if cost_rows[node_signature(root)][_ROOT_STATE] == _INFINITE_COST:
        raise ParsimonyAnalysisError(
            "camin-sokal scoring could not realize the observed binary pattern without reversals",
            code="parsimony_irreversible_pattern_impossible",
            details={"character_id": character_id},
        )
    branch_changes: list[CaminSokalBranchChange] = []

    def assign(node: TreeNode, parent_state: str) -> None:
        for child in node.children:
            child_state = _choose_child_state(
                parent_state,
                child,
                cost_rows=cost_rows,
            )
            if parent_state == "0" and child_state == "1":
                branch_changes.append(
                    CaminSokalBranchChange(
                        character_id=character_id,
                        change_kind="gain",
                        node=node_signature(child),
                        node_name=child.name,
                        descendant_taxa=node_descendant_taxa(child),
                    )
                )
            assign(child, child_state)

    assign(root, _ROOT_STATE)
    return branch_changes


def _compute_irreversible_costs(
    root: TreeNode,
    *,
    character_id: str,
    states_by_taxon: dict[str, dict[str, str]],
) -> dict[str, dict[str, float]]:
    cost_rows: dict[str, dict[str, float]] = {}
    for node in root.iter_nodes(order="postorder"):
        node_key = node_signature(node)
        if node.is_leaf():
            if node.name is None:
                raise ParsimonyAnalysisError(
                    "camin-sokal scoring requires every leaf to have a taxon label",
                    code="parsimony_tree_unnamed_tip",
                )
            observed_state = states_by_taxon[node.name][character_id]
            cost_rows[node_key] = {
                state: 0.0 if state == observed_state else _INFINITE_COST
                for state in _STATES
            }
            continue
        node_costs: dict[str, float] = {}
        for state in _STATES:
            total_cost = 0.0
            for child in node.children:
                child_key = node_signature(child)
                candidate_costs = [
                    cost_rows[child_key][child_state]
                    + _transition_cost(state, child_state)
                    for child_state in _STATES
                ]
                total_cost += min(candidate_costs)
            node_costs[state] = total_cost
        cost_rows[node_key] = node_costs
    return cost_rows


def _choose_child_state(
    parent_state: str,
    child: TreeNode,
    *,
    cost_rows: dict[str, dict[str, float]],
) -> str:
    child_key = node_signature(child)
    best_cost = min(
        cost_rows[child_key][child_state] + _transition_cost(parent_state, child_state)
        for child_state in _STATES
    )
    preferred_state_order = (parent_state, "0", "1")
    for child_state in preferred_state_order:
        total_cost = cost_rows[child_key][child_state] + _transition_cost(
            parent_state,
            child_state,
        )
        if total_cost == best_cost:
            return child_state
    raise AssertionError("failed to select one irreversible child state")


def _transition_cost(parent_state: str, child_state: str) -> float:
    if parent_state == child_state:
        return 0.0
    if parent_state == "0" and child_state == "1":
        return 1.0
    return _INFINITE_COST
