from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .fitch import _leaf_taxa, _resolve_tree
from .matrix import load_parsimony_character_matrix
from .models import (
    DolloBranchChange,
    DolloCharacterScore,
    DolloScoreReport,
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
)
from .weights import resolve_parsimony_character_weights

_BINARY_STATES = frozenset({"0", "1"})


def score_dollo(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
) -> DolloScoreReport:
    """Score one binary character matrix on one tree with Dollo parsimony."""
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
            "dollo scoring requires every tree taxon to be present in the character matrix",
            code="parsimony_matrix_missing_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "missing_taxa": missing_from_matrix,
            },
        )
    if extra_in_matrix:
        raise ParsimonyAnalysisError(
            "dollo scoring requires matrix taxa to match the tree tips exactly",
            code="parsimony_matrix_extra_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "extra_taxa": extra_in_matrix,
            },
        )

    step_rows: list[DolloCharacterScore] = []
    branch_change_rows: list[DolloBranchChange] = []
    total_gains = 0
    total_losses = 0
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
                "dollo scoring requires binary 0/1 character states unless the trait is binarized first",
                code="parsimony_matrix_multistate_not_binarized",
                details={
                    "character_id": character_id,
                    "invalid_states": invalid_states,
                },
            )
        derived_taxa = sorted(
            taxon
            for taxon, states_by_character in resolved_matrix.states_by_taxon.items()
            if states_by_character[character_id] == "1"
        )
        if not derived_taxa:
            step_rows.append(
                DolloCharacterScore(
                    character_id=character_id,
                    step_count=0,
                    derived_taxon_count=0,
                    gain_node=None,
                    gain_node_name=None,
                    gain_descendant_taxa=[],
                    total_losses=0,
                    impossible_state_warning=None,
                    character_weight=character_weight,
                    weighted_score=0.0,
                )
            )
            continue

        gain_node = _lowest_common_ancestor_for_taxa(resolved_tree, derived_taxa)
        total_gains += 1
        gain_signature = node_signature(gain_node)
        gain_descendant_taxa = node_descendant_taxa(gain_node)
        branch_change_rows.append(
            DolloBranchChange(
                character_id=character_id,
                change_kind="gain",
                node=gain_signature,
                node_name=gain_node.name,
                descendant_taxa=gain_descendant_taxa,
            )
        )
        loss_nodes = _loss_nodes_for_gain(gain_node, derived_taxa=derived_taxa)
        total_losses += len(loss_nodes)
        step_count = 1 + len(loss_nodes)
        weighted_score = step_count * character_weight
        total_weighted_score += weighted_score
        for loss_node in loss_nodes:
            branch_change_rows.append(
                DolloBranchChange(
                    character_id=character_id,
                    change_kind="loss",
                    node=node_signature(loss_node),
                    node_name=loss_node.name,
                    descendant_taxa=node_descendant_taxa(loss_node),
                )
            )
        step_rows.append(
            DolloCharacterScore(
                character_id=character_id,
                step_count=step_count,
                derived_taxon_count=len(derived_taxa),
                gain_node=gain_signature,
                gain_node_name=gain_node.name,
                gain_descendant_taxa=gain_descendant_taxa,
                total_losses=len(loss_nodes),
                impossible_state_warning=None,
                character_weight=character_weight,
                weighted_score=weighted_score,
            )
        )
    return DolloScoreReport(
        algorithm="dollo",
        tree_path=tree_path,
        matrix_path=matrix_path,
        taxon_column=resolved_matrix.taxon_column,
        taxon_count=len(leaf_taxa),
        character_count=resolved_matrix.character_count,
        total_gains=total_gains,
        total_losses=total_losses,
        weights_path=resolved_weights.weights_path,
        total_weighted_score=total_weighted_score,
        step_rows=step_rows,
        branch_change_rows=branch_change_rows,
    )


def _lowest_common_ancestor_for_taxa(
    tree: PhyloTree,
    taxa: list[str],
) -> TreeNode:
    target_taxa = set(taxa)
    best_node = tree.root
    best_size = len(node_descendant_taxa(tree.root))
    for node in tree.iter_nodes(order="preorder"):
        descendants = set(node_descendant_taxa(node))
        if target_taxa.issubset(descendants) and len(descendants) <= best_size:
            best_node = node
            best_size = len(descendants)
    return best_node


def _loss_nodes_for_gain(
    gain_node: TreeNode, *, derived_taxa: list[str]
) -> list[TreeNode]:
    target_taxa = set(derived_taxa)
    loss_nodes: list[TreeNode] = []

    def visit(node: TreeNode) -> None:
        for child in node.children:
            descendants = set(node_descendant_taxa(child))
            if descendants.isdisjoint(target_taxa):
                loss_nodes.append(child)
                continue
            visit(child)

    visit(gain_node)
    return loss_nodes
