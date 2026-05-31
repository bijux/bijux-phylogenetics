from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .fitch import _leaf_taxa, _resolve_tree
from .matrix import load_parsimony_character_matrix
from .models import (
    FitchCharacterMatrix,
    ParsimonyAncestralState,
    ParsimonyCharacterWeights,
    ParsimonyReconstructionBranchChange,
    ParsimonyReconstructionCharacterScore,
    ParsimonyReconstructionNodeState,
    ParsimonyReconstructionReport,
)
from .weights import resolve_parsimony_character_weights


def reconstruct_acctran(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
) -> ParsimonyReconstructionReport:
    """Resolve unordered parsimony ambiguities toward earlier branches with ACCTRAN."""
    return _reconstruct_parsimony_report(
        tree,
        matrix,
        algorithm="acctran",
        character_weights=character_weights,
    )


def reconstruct_deltran(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
) -> ParsimonyReconstructionReport:
    """Resolve unordered parsimony ambiguities toward later branches with DELTRAN."""
    return _reconstruct_parsimony_report(
        tree,
        matrix,
        algorithm="deltran",
        character_weights=character_weights,
    )


def _reconstruct_deltran(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
) -> ParsimonyReconstructionReport:
    """Backward-compatible internal alias for DELTRAN reconstruction."""
    return reconstruct_deltran(tree, matrix, character_weights=character_weights)


def _reconstruct_parsimony_report(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    algorithm: str,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
) -> ParsimonyReconstructionReport:
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
            f"{algorithm} reconstruction requires every tree taxon to be present in the character matrix",
            code="parsimony_matrix_missing_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "missing_taxa": missing_from_matrix,
            },
        )
    if extra_in_matrix:
        raise ParsimonyAnalysisError(
            f"{algorithm} reconstruction requires matrix taxa to match the tree tips exactly",
            code="parsimony_matrix_extra_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "extra_taxa": extra_in_matrix,
            },
        )

    node_depths = _node_depths(resolved_tree.root)
    step_rows: list[ParsimonyReconstructionCharacterScore] = []
    node_state_rows: list[ParsimonyReconstructionNodeState] = []
    ancestral_state_rows: list[ParsimonyAncestralState] = []
    branch_change_rows: list[ParsimonyReconstructionBranchChange] = []
    total_steps = 0
    total_weighted_score = 0.0
    preorder_nodes = list(resolved_tree.iter_nodes(order="preorder"))
    for character_id in resolved_matrix.character_ids:
        observed_states = sorted(
            {
                states_by_character[character_id]
                for states_by_character in resolved_matrix.states_by_taxon.values()
            }
        )
        dynamic_rows = _compute_dynamic_rows(
            resolved_tree.root,
            character_id=character_id,
            states_by_taxon=resolved_matrix.states_by_taxon,
            states=observed_states,
            node_depths=node_depths,
            algorithm=algorithm,
        )
        possible_state_sets = _possible_state_sets(
            resolved_tree.root,
            character_id=character_id,
            states_by_taxon=resolved_matrix.states_by_taxon,
        )
        root_signature = node_signature(resolved_tree.root)
        root_state = min(
            observed_states,
            key=lambda state: (
                dynamic_rows[root_signature][state][0],
                dynamic_rows[root_signature][state][1],
                state,
            ),
        )
        step_count = int(dynamic_rows[root_signature][root_state][0])
        total_steps += step_count
        character_weight = resolved_weights.weights_by_character[character_id]
        weighted_score = step_count * character_weight
        total_weighted_score += weighted_score
        assigned_states = {root_signature: root_state}
        _assign_child_states(
            resolved_tree.root,
            parent_state=root_state,
            character_id=character_id,
            states=observed_states,
            dynamic_rows=dynamic_rows,
            node_depths=node_depths,
            algorithm=algorithm,
            assigned_states=assigned_states,
            branch_change_rows=branch_change_rows,
        )
        step_rows.append(
            ParsimonyReconstructionCharacterScore(
                character_id=character_id,
                step_count=step_count,
                observed_states=observed_states,
                root_state=root_state,
                character_weight=character_weight,
                weighted_score=weighted_score,
            )
        )
        for node in preorder_nodes:
            signature = node_signature(node)
            descendant_taxa = node_descendant_taxa(node)
            node_state_rows.append(
                ParsimonyReconstructionNodeState(
                    character_id=character_id,
                    node=signature,
                    node_name=node.name,
                    descendant_taxa=descendant_taxa,
                    resolved_state=assigned_states[signature],
                    is_tip=node.is_leaf(),
                    observed_state=(
                        resolved_matrix.states_by_taxon[node.name][character_id]
                        if node.is_leaf() and node.name is not None
                        else None
                    ),
                )
            )
            if not node.is_leaf():
                ancestral_state_rows.append(
                    ParsimonyAncestralState(
                        node_id=signature,
                        clade_id="|".join(descendant_taxa),
                        character_id=character_id,
                        possible_states=sorted(possible_state_sets[signature]),
                        chosen_state=assigned_states[signature],
                        method=algorithm,
                        ambiguous=len(possible_state_sets[signature]) > 1,
                    )
                )
    return ParsimonyReconstructionReport(
        algorithm=algorithm,
        tree_path=tree_path,
        matrix_path=matrix_path,
        taxon_column=resolved_matrix.taxon_column,
        taxon_count=len(leaf_taxa),
        character_count=resolved_matrix.character_count,
        total_steps=total_steps,
        weights_path=resolved_weights.weights_path,
        total_weighted_score=total_weighted_score,
        step_rows=step_rows,
        node_state_rows=node_state_rows,
        ancestral_state_rows=ancestral_state_rows,
        branch_change_rows=branch_change_rows,
    )


def _compute_dynamic_rows(
    root: TreeNode,
    *,
    character_id: str,
    states_by_taxon: dict[str, dict[str, str]],
    states: list[str],
    node_depths: dict[str, int],
    algorithm: str,
) -> dict[str, dict[str, tuple[int, int]]]:
    dynamic_rows: dict[str, dict[str, tuple[int, int]]] = {}
    for node in root.iter_nodes(order="postorder"):
        signature = node_signature(node)
        if node.is_leaf():
            if node.name is None:
                raise ParsimonyAnalysisError(
                    f"{algorithm} reconstruction requires every leaf to have a taxon label",
                    code="parsimony_tree_unnamed_tip",
                )
            observed_state = states_by_taxon[node.name][character_id]
            dynamic_rows[signature] = {
                state: (0, 0) if state == observed_state else (10**9, 10**9)
                for state in states
            }
            continue
        state_rows: dict[str, tuple[int, int]] = {}
        for state in states:
            total_steps = 0
            total_secondary = 0
            for child in node.children:
                child_signature = node_signature(child)
                best_child_state = min(
                    states,
                    key=lambda child_state: _child_objective(
                        parent_state=state,
                        child_state=child_state,
                        child_signature=child_signature,
                        dynamic_rows=dynamic_rows,
                        node_depths=node_depths,
                        algorithm=algorithm,
                    ),
                )
                child_steps, child_secondary, _ = _child_objective(
                    parent_state=state,
                    child_state=best_child_state,
                    child_signature=child_signature,
                    dynamic_rows=dynamic_rows,
                    node_depths=node_depths,
                    algorithm=algorithm,
                )
                total_steps += child_steps
                total_secondary += child_secondary
            state_rows[state] = (total_steps, total_secondary)
        dynamic_rows[signature] = state_rows
    return dynamic_rows


def _possible_state_sets(
    root: TreeNode,
    *,
    character_id: str,
    states_by_taxon: dict[str, dict[str, str]],
) -> dict[str, set[str]]:
    candidate_sets: dict[str, set[str]] = {}
    for node in root.iter_nodes(order="postorder"):
        signature = node_signature(node)
        if node.is_leaf():
            if node.name is None:
                raise ParsimonyAnalysisError(
                    "parsimony reconstruction requires every leaf to have a taxon label",
                    code="parsimony_tree_unnamed_tip",
                )
            candidate_sets[signature] = {states_by_taxon[node.name][character_id]}
            continue
        child_sets = [
            set(candidate_sets[node_signature(child)]) for child in node.children
        ]
        candidate = child_sets[0]
        for child_set in child_sets[1:]:
            intersection = candidate & child_set
            if intersection:
                candidate = intersection
            else:
                candidate |= child_set
        candidate_sets[signature] = candidate
    return candidate_sets


def _assign_child_states(
    node: TreeNode,
    *,
    parent_state: str,
    character_id: str,
    states: list[str],
    dynamic_rows: dict[str, dict[str, tuple[int, int]]],
    node_depths: dict[str, int],
    algorithm: str,
    assigned_states: dict[str, str],
    branch_change_rows: list[ParsimonyReconstructionBranchChange],
) -> None:
    for child in node.children:
        child_signature = node_signature(child)
        possible_states = _possible_child_states(
            parent_state=parent_state,
            child_signature=child_signature,
            states=states,
            dynamic_rows=dynamic_rows,
            node_depths=node_depths,
            algorithm=algorithm,
        )
        child_state = min(
            possible_states,
            key=lambda state: _child_objective(
                parent_state=parent_state,
                child_state=state,
                child_signature=child_signature,
                dynamic_rows=dynamic_rows,
                node_depths=node_depths,
                algorithm=algorithm,
            ),
        )
        assigned_states[child_signature] = child_state
        if child_state != parent_state:
            branch_change_rows.append(
                ParsimonyReconstructionBranchChange(
                    character_id=character_id,
                    branch_id=child_signature,
                    parent_node=node_signature(node),
                    parent_state=parent_state,
                    child_node=child_signature,
                    child_node_name=child.name,
                    child_descendant_taxa=node_descendant_taxa(child),
                    change_from=parent_state,
                    change_to=child_state,
                    ambiguous=len(possible_states) > 1,
                )
            )
        _assign_child_states(
            child,
            parent_state=child_state,
            character_id=character_id,
            states=states,
            dynamic_rows=dynamic_rows,
            node_depths=node_depths,
            algorithm=algorithm,
            assigned_states=assigned_states,
            branch_change_rows=branch_change_rows,
        )


def _child_objective(
    *,
    parent_state: str,
    child_state: str,
    child_signature: str,
    dynamic_rows: dict[str, dict[str, tuple[int, int]]],
    node_depths: dict[str, int],
    algorithm: str,
) -> tuple[int, int, str]:
    child_steps, child_secondary = dynamic_rows[child_signature][child_state]
    branch_change = 0 if parent_state == child_state else 1
    branch_depth = node_depths[child_signature]
    branch_secondary = branch_depth if algorithm == "acctran" else -branch_depth
    return (
        child_steps + branch_change,
        child_secondary + (branch_secondary if branch_change else 0),
        child_state,
    )


def _possible_child_states(
    *,
    parent_state: str,
    child_signature: str,
    states: list[str],
    dynamic_rows: dict[str, dict[str, tuple[int, int]]],
    node_depths: dict[str, int],
    algorithm: str,
) -> list[str]:
    best_objective = min(
        _child_objective(
            parent_state=parent_state,
            child_state=state,
            child_signature=child_signature,
            dynamic_rows=dynamic_rows,
            node_depths=node_depths,
            algorithm=algorithm,
        )[:2]
        for state in states
    )
    return [
        state
        for state in states
        if _child_objective(
            parent_state=parent_state,
            child_state=state,
            child_signature=child_signature,
            dynamic_rows=dynamic_rows,
            node_depths=node_depths,
            algorithm=algorithm,
        )[:2]
        == best_objective
    ]


def _node_depths(root: TreeNode) -> dict[str, int]:
    depths: dict[str, int] = {}

    def visit(node: TreeNode, depth: int) -> None:
        depths[node_signature(node)] = depth
        for child in node.children:
            visit(child, depth + 1)

    visit(root, 0)
    return depths
