from __future__ import annotations

from collections.abc import Mapping
import math
from pathlib import Path

from bijux_phylogenetics.io.newick import load_newick
from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode, descendant_taxa
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

from .cost_matrix import load_sankoff_cost_matrix
from .exact_resampling import SUPPORTED_RESAMPLED_PARSIMONY_METHODS
from .matrix import load_parsimony_character_matrix
from .models import FitchCharacterMatrix, ParsimonyCharacterWeights, SankoffCostMatrix
from .tree_length import tree_length
from .weights import resolve_parsimony_character_weights


def resolve_topology_search_tree(
    tree: PhyloTree | Path,
) -> tuple[PhyloTree, Path | None]:
    """Resolve one search tree input and preserve its file path when present."""
    resolved_tree_path = tree if isinstance(tree, Path) else None
    resolved_tree = load_newick(tree) if isinstance(tree, Path) else tree.copy()
    return resolved_tree.refresh(), resolved_tree_path


def resolve_topology_search_matrix(
    matrix: FitchCharacterMatrix | Path,
    *,
    taxon_column: str | None,
) -> FitchCharacterMatrix:
    """Resolve one topology-search character matrix from disk or memory."""
    return (
        matrix
        if isinstance(matrix, FitchCharacterMatrix)
        else load_parsimony_character_matrix(matrix, taxon_column=taxon_column)
    )


def resolve_topology_search_method(method: str, *, workflow_name: str) -> str:
    """Validate one supported parsimony method for topology search."""
    resolved_method = method.strip().lower()
    if resolved_method not in SUPPORTED_RESAMPLED_PARSIMONY_METHODS:
        raise ParsimonyAnalysisError(
            f"{workflow_name} requires one supported parsimony method",
            code=f"{_workflow_slug(workflow_name)}_method_unsupported",
            details={
                "method": method,
                "supported_methods": sorted(SUPPORTED_RESAMPLED_PARSIMONY_METHODS),
            },
        )
    return resolved_method


def resolve_topology_search_cost_matrix(
    *,
    method: str,
    cost_matrix: SankoffCostMatrix | Path | None,
    allow_asymmetric_costs: bool,
    workflow_name: str,
) -> SankoffCostMatrix | None:
    """Resolve one optional Sankoff cost matrix for topology search."""
    if method != "sankoff":
        return None
    if cost_matrix is None:
        raise ParsimonyAnalysisError(
            f"{workflow_name} requires a Sankoff cost matrix when method='sankoff'",
            code=f"{_workflow_slug(workflow_name)}_cost_matrix_required",
            details={"method": method},
        )
    if isinstance(cost_matrix, SankoffCostMatrix):
        return cost_matrix
    return load_sankoff_cost_matrix(
        cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
    )


def resolve_topology_search_weights(
    matrix: FitchCharacterMatrix,
    character_weights: ParsimonyCharacterWeights | Mapping[str, float] | Path | None,
) -> ParsimonyCharacterWeights:
    """Resolve one optional explicit weight table for topology search."""
    return resolve_parsimony_character_weights(
        matrix.character_ids,
        character_weights,
    )


def validate_topology_search_tree(tree: PhyloTree, *, workflow_name: str) -> None:
    """Require a structurally valid strictly binary rooted tree."""
    errors = tree.validation_errors()
    if errors:
        raise ParsimonyAnalysisError(
            f"{workflow_name} requires a structurally valid tree",
            code=f"{_workflow_slug(workflow_name)}_tree_invalid",
            details={"validation_errors": errors},
        )
    if len(tree.root.children) != 2:
        raise ParsimonyAnalysisError(
            f"{workflow_name} requires a rooted binary tree",
            code=f"{_workflow_slug(workflow_name)}_root_not_binary",
            details={"root_child_count": len(tree.root.children)},
        )
    invalid_internal_nodes = [
        node.node_id
        for node in tree.iter_internal_nodes(order="preorder")
        if len(node.children) != 2
    ]
    if invalid_internal_nodes:
        raise ParsimonyAnalysisError(
            f"{workflow_name} requires a strictly binary tree",
            code=f"{_workflow_slug(workflow_name)}_tree_not_binary",
            details={"invalid_internal_nodes": invalid_internal_nodes},
        )


def score_topology_search_tree(
    tree: PhyloTree,
    matrix: FitchCharacterMatrix,
    *,
    method: str,
    state_order: list[str] | None,
    cost_matrix: SankoffCostMatrix | None,
    allow_asymmetric_costs: bool,
    character_weights: ParsimonyCharacterWeights,
) -> float:
    """Score one topology-search tree through the governed parsimony scorer surface."""
    return tree_length(
        tree,
        matrix,
        method=method,
        state_order=state_order,
        cost_matrix=cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        character_weights=character_weights,
    ).total_score


def topology_search_prefer_score(
    left_score: float,
    left_newick: str,
    right_score: float,
    right_newick: str,
) -> bool:
    """Prefer lower scores, then deterministic Newick order, across search moves."""
    if left_score < right_score and not math.isclose(left_score, right_score):
        return True
    if right_score < left_score and not math.isclose(left_score, right_score):
        return False
    return left_newick < right_newick


def topology_search_node_sort_key(node: TreeNode) -> tuple[int, tuple[str, ...]]:
    """Sort candidate search branches deterministically by descendant taxa."""
    descendants = tuple(descendant_taxa(node))
    return (len(descendants), descendants)


def topology_search_clade_id(node: TreeNode) -> str:
    """Render one branch endpoint as a stable descendant-clade identifier."""
    return canonical_clade_id(frozenset(descendant_taxa(node)))


def require_topology_search_node_id(node: TreeNode) -> str:
    """Require one refreshed stable node identifier on a search tree."""
    if node.node_id is None:
        raise AssertionError("topology search requires refreshed node identities")
    return node.node_id


def _workflow_slug(workflow_name: str) -> str:
    return workflow_name.strip().lower().replace(" ", "_")
