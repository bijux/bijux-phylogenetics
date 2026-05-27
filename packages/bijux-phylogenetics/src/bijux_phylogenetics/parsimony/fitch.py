from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .matrix import load_fitch_character_matrix
from .models import (
    FitchCharacterMatrix,
    FitchCharacterScore,
    FitchNodeStateSet,
    FitchScoreReport,
    ParsimonyCharacterWeights,
)
from .weights import resolve_parsimony_character_weights


def score_fitch(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
) -> FitchScoreReport:
    """Score one taxon-by-character matrix on one tree with unordered Fitch parsimony."""
    resolved_tree, tree_path = _resolve_tree(tree)
    resolved_matrix = (
        matrix
        if isinstance(matrix, FitchCharacterMatrix)
        else load_fitch_character_matrix(matrix)
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
            "unordered fitch scoring requires every tree taxon to be present in the character matrix",
            code="parsimony_matrix_missing_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "missing_taxa": missing_from_matrix,
            },
        )
    if extra_in_matrix:
        raise ParsimonyAnalysisError(
            "unordered fitch scoring requires matrix taxa to match the tree tips exactly",
            code="parsimony_matrix_extra_taxa",
            details={
                "tree_path": None if tree_path is None else str(tree_path),
                "matrix_path": None if matrix_path is None else str(matrix_path),
                "extra_taxa": extra_in_matrix,
            },
        )

    step_rows: list[FitchCharacterScore] = []
    node_state_rows: list[FitchNodeStateSet] = []
    total_steps = 0
    total_weighted_score = 0.0
    preorder_nodes = list(resolved_tree.iter_nodes(order="preorder"))
    for character_id in resolved_matrix.character_ids:
        candidate_sets, step_count = _score_character(
            resolved_tree.root,
            character_id=character_id,
            states_by_taxon=resolved_matrix.states_by_taxon,
        )
        total_steps += step_count
        character_weight = resolved_weights.weights_by_character[character_id]
        weighted_score = step_count * character_weight
        total_weighted_score += weighted_score
        observed_states = sorted(
            {
                states_by_character[character_id]
                for states_by_character in resolved_matrix.states_by_taxon.values()
            }
        )
        step_rows.append(
            FitchCharacterScore(
                character_id=character_id,
                step_count=step_count,
                observed_states=observed_states,
                character_weight=character_weight,
                weighted_score=weighted_score,
            )
        )
        for node in preorder_nodes:
            node_key = node.node_id or node_signature(node)
            node_state_rows.append(
                FitchNodeStateSet(
                    character_id=character_id,
                    node=node_signature(node),
                    node_name=node.name,
                    descendant_taxa=node_descendant_taxa(node),
                    state_set=sorted(candidate_sets[node_key]),
                    is_tip=node.is_leaf(),
                    observed_state=(
                        resolved_matrix.states_by_taxon[node.name][character_id]
                        if node.is_leaf() and node.name is not None
                        else None
                    ),
                )
            )
    return FitchScoreReport(
        algorithm="unordered-fitch",
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
    )


def _resolve_tree(tree: PhyloTree | Path) -> tuple[PhyloTree, Path | None]:
    if isinstance(tree, Path):
        return load_tree(tree), tree
    return tree, None


def _leaf_taxa(tree: PhyloTree) -> list[str]:
    taxa = [node.name for node in tree.iter_leaves()]
    unnamed_taxa = [taxon for taxon in taxa if taxon is None]
    if unnamed_taxa:
        raise ParsimonyAnalysisError(
            "unordered fitch scoring requires every tree tip to have a taxon label",
            code="parsimony_tree_unnamed_tip",
        )
    return sorted(taxon for taxon in taxa if taxon is not None)


def _score_character(
    root: TreeNode,
    *,
    character_id: str,
    states_by_taxon: dict[str, dict[str, str]],
) -> tuple[dict[str, set[str]], int]:
    candidate_sets: dict[str, set[str]] = {}
    step_count = 0
    for node in root.iter_nodes(order="postorder"):
        node_key = node.node_id or node_signature(node)
        if node.is_leaf():
            if node.name is None:
                raise ParsimonyAnalysisError(
                    "unordered fitch scoring requires every leaf to have a taxon label",
                    code="parsimony_tree_unnamed_tip",
                )
            candidate_sets[node_key] = {states_by_taxon[node.name][character_id]}
            continue
        child_sets = [
            set(candidate_sets[child.node_id or node_signature(child)])
            for child in node.children
        ]
        candidate = child_sets[0]
        for child_set in child_sets[1:]:
            intersection = candidate & child_set
            if intersection:
                candidate = intersection
            else:
                candidate |= child_set
                step_count += 1
        candidate_sets[node_key] = candidate
    return candidate_sets, step_count
