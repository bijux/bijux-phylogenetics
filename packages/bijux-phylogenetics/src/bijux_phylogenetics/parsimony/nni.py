from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick, load_newick
from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode, descendant_taxa

from .cost_matrix import load_sankoff_cost_matrix
from .exact_resampling import SUPPORTED_RESAMPLED_PARSIMONY_METHODS
from .matrix import load_parsimony_character_matrix
from .models import (
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
    ParsimonyNniSearchReport,
    ParsimonyNniTraceRow,
    SankoffCostMatrix,
)
from .tree_length import tree_length
from .weights import resolve_parsimony_character_weights
from ..runtime.errors import ParsimonyAnalysisError


@dataclass(frozen=True, slots=True)
class _NniMoveCandidate:
    parent_node_id: str
    child_node_id: str
    sibling_node_id: str
    exchanged_child_node_id: str
    pivot_branch_id: str
    sibling_clade_id: str
    exchanged_clade_id: str


@dataclass(frozen=True, slots=True)
class _ScoredNniMove:
    candidate: _NniMoveCandidate
    tree: PhyloTree
    score: float
    newick: str


def search_parsimony_nni(
    tree: PhyloTree | Path,
    matrix: FitchCharacterMatrix | Path,
    *,
    method: str,
    taxon_column: str | None = None,
    state_order: list[str] | None = None,
    cost_matrix: SankoffCostMatrix | Path | None = None,
    allow_asymmetric_costs: bool = False,
    character_weights: (
        ParsimonyCharacterWeights | Mapping[str, float] | Path | None
    ) = None,
) -> ParsimonyNniSearchReport:
    """Search one rooted binary tree by accepting score-improving rooted NNI moves."""
    resolved_tree_path = tree if isinstance(tree, Path) else None
    resolved_tree = _resolve_tree(tree)
    _validate_nni_tree(resolved_tree)
    resolved_matrix = (
        matrix
        if isinstance(matrix, FitchCharacterMatrix)
        else load_parsimony_character_matrix(matrix, taxon_column=taxon_column)
    )
    resolved_method = _resolve_method(method)
    resolved_cost_matrix = _resolve_cost_matrix(
        method=resolved_method,
        cost_matrix=cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
    )
    resolved_weights = resolve_parsimony_character_weights(
        resolved_matrix.character_ids,
        character_weights,
    )
    current_tree = resolved_tree.copy().refresh()
    current_score = _score_tree(
        current_tree,
        resolved_matrix,
        method=resolved_method,
        state_order=state_order,
        cost_matrix=resolved_cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        character_weights=resolved_weights,
    )
    start_tree_newick = dumps_newick(current_tree)
    trace_rows = [
        ParsimonyNniTraceRow(
            event_index=1,
            event_kind="start",
            iteration=0,
            score_before=None,
            score_after=current_score,
            score_delta=None,
            tree_before_newick=None,
            tree_after_newick=start_tree_newick,
            pivot_branch_id=None,
            sibling_clade_id=None,
            exchanged_clade_id=None,
            stopping_reason=None,
        )
    ]
    accepted_move_count = 0
    evaluated_neighbor_count = 0
    while True:
        improving_move: _ScoredNniMove | None = None
        for candidate in _iter_nni_move_candidates(current_tree):
            neighbor_tree = _apply_nni_move(current_tree, candidate)
            neighbor_score = _score_tree(
                neighbor_tree,
                resolved_matrix,
                method=resolved_method,
                state_order=state_order,
                cost_matrix=resolved_cost_matrix,
                allow_asymmetric_costs=allow_asymmetric_costs,
                character_weights=resolved_weights,
            )
            evaluated_neighbor_count += 1
            neighbor_newick = dumps_newick(neighbor_tree)
            if neighbor_score > current_score or math.isclose(
                neighbor_score,
                current_score,
            ):
                continue
            scored_move = _ScoredNniMove(
                candidate=candidate,
                tree=neighbor_tree,
                score=neighbor_score,
                newick=neighbor_newick,
            )
            if improving_move is None or _prefer_scored_move(
                scored_move,
                improving_move,
            ):
                improving_move = scored_move
        if improving_move is None:
            stopping_reason = "no-improving-neighbor"
            break
        accepted_move_count += 1
        tree_before_newick = dumps_newick(current_tree)
        score_before = current_score
        current_tree = improving_move.tree.copy().refresh()
        current_score = improving_move.score
        trace_rows.append(
            ParsimonyNniTraceRow(
                event_index=len(trace_rows) + 1,
                event_kind="accepted-move",
                iteration=accepted_move_count,
                score_before=score_before,
                score_after=current_score,
                score_delta=current_score - score_before,
                tree_before_newick=tree_before_newick,
                tree_after_newick=improving_move.newick,
                pivot_branch_id=improving_move.candidate.pivot_branch_id,
                sibling_clade_id=improving_move.candidate.sibling_clade_id,
                exchanged_clade_id=improving_move.candidate.exchanged_clade_id,
                stopping_reason=None,
            )
        )
    trace_rows.append(
        ParsimonyNniTraceRow(
            event_index=len(trace_rows) + 1,
            event_kind="final",
            iteration=accepted_move_count,
            score_before=None,
            score_after=current_score,
            score_delta=None,
            tree_before_newick=None,
            tree_after_newick=dumps_newick(current_tree),
            pivot_branch_id=None,
            sibling_clade_id=None,
            exchanged_clade_id=None,
            stopping_reason=stopping_reason,
        )
    )
    return ParsimonyNniSearchReport(
        algorithm="parsimony-nni-search",
        method=resolved_method,
        tree_path=resolved_tree_path,
        matrix_path=resolved_matrix.matrix_path,
        cost_matrix_path=None
        if resolved_cost_matrix is None
        else resolved_cost_matrix.matrix_path,
        weights_path=resolved_weights.weights_path,
        taxon_column=resolved_matrix.taxon_column,
        taxon_count=resolved_matrix.taxon_count,
        character_count=resolved_matrix.character_count,
        start_tree_newick=start_tree_newick,
        start_score=trace_rows[0].score_after,
        final_tree_newick=trace_rows[-1].tree_after_newick,
        final_score=trace_rows[-1].score_after,
        accepted_move_count=accepted_move_count,
        evaluated_neighbor_count=evaluated_neighbor_count,
        stopping_reason=stopping_reason,
        trace_rows=trace_rows,
    )


def _resolve_tree(tree: PhyloTree | Path) -> PhyloTree:
    resolved_tree = load_newick(tree) if isinstance(tree, Path) else tree.copy()
    return resolved_tree.refresh()


def _resolve_method(method: str) -> str:
    resolved_method = method.strip().lower()
    if resolved_method not in SUPPORTED_RESAMPLED_PARSIMONY_METHODS:
        raise ParsimonyAnalysisError(
            "parsimony NNI search requires one supported parsimony method",
            code="parsimony_nni_method_unsupported",
            details={
                "method": method,
                "supported_methods": sorted(SUPPORTED_RESAMPLED_PARSIMONY_METHODS),
            },
        )
    return resolved_method


def _resolve_cost_matrix(
    *,
    method: str,
    cost_matrix: SankoffCostMatrix | Path | None,
    allow_asymmetric_costs: bool,
) -> SankoffCostMatrix | None:
    if method != "sankoff":
        return None
    if cost_matrix is None:
        raise ParsimonyAnalysisError(
            "parsimony NNI search requires a Sankoff cost matrix when method='sankoff'",
            code="parsimony_nni_cost_matrix_required",
            details={"method": method},
        )
    if isinstance(cost_matrix, SankoffCostMatrix):
        return cost_matrix
    return load_sankoff_cost_matrix(
        cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
    )


def _validate_nni_tree(tree: PhyloTree) -> None:
    errors = tree.validation_errors()
    if errors:
        raise ParsimonyAnalysisError(
            "parsimony NNI search requires a structurally valid tree",
            code="parsimony_nni_tree_invalid",
            details={"validation_errors": errors},
        )
    if len(tree.root.children) != 2:
        raise ParsimonyAnalysisError(
            "parsimony NNI search requires a rooted binary tree",
            code="parsimony_nni_root_not_binary",
            details={"root_child_count": len(tree.root.children)},
        )
    invalid_internal_nodes = [
        node.node_id
        for node in tree.iter_internal_nodes(order="preorder")
        if len(node.children) != 2
    ]
    if invalid_internal_nodes:
        raise ParsimonyAnalysisError(
            "parsimony NNI search requires a strictly binary tree",
            code="parsimony_nni_tree_not_binary",
            details={"invalid_internal_nodes": invalid_internal_nodes},
        )


def _score_tree(
    tree: PhyloTree,
    matrix: FitchCharacterMatrix,
    *,
    method: str,
    state_order: list[str] | None,
    cost_matrix: SankoffCostMatrix | None,
    allow_asymmetric_costs: bool,
    character_weights: ParsimonyCharacterWeights,
) -> float:
    return tree_length(
        tree,
        matrix,
        method=method,
        state_order=state_order,
        cost_matrix=cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        character_weights=character_weights,
    ).total_score


def _iter_nni_move_candidates(tree: PhyloTree):
    for parent in tree.iter_internal_nodes(order="preorder"):
        if len(parent.children) != 2:
            continue
        sorted_parent_children = sorted(parent.children, key=_node_sort_key)
        for child in sorted_parent_children:
            if child.is_leaf() or len(child.children) != 2:
                continue
            sibling = next(
                candidate
                for candidate in sorted_parent_children
                if candidate is not child
            )
            sorted_child_children = sorted(child.children, key=_node_sort_key)
            for exchanged_child in sorted_child_children:
                yield _NniMoveCandidate(
                    parent_node_id=_require_node_id(parent),
                    child_node_id=_require_node_id(child),
                    sibling_node_id=_require_node_id(sibling),
                    exchanged_child_node_id=_require_node_id(exchanged_child),
                    pivot_branch_id=_clade_id(child),
                    sibling_clade_id=_clade_id(sibling),
                    exchanged_clade_id=_clade_id(exchanged_child),
                )


def _apply_nni_move(tree: PhyloTree, candidate: _NniMoveCandidate) -> PhyloTree:
    swapped_tree = tree.copy().refresh()
    parent = swapped_tree.node_by_id(candidate.parent_node_id)
    child = swapped_tree.node_by_id(candidate.child_node_id)
    sibling = swapped_tree.node_by_id(candidate.sibling_node_id)
    exchanged_child = swapped_tree.node_by_id(candidate.exchanged_child_node_id)
    remaining_child = next(
        branch for branch in child.children if branch is not exchanged_child
    )
    child.replace_children([remaining_child, sibling])
    parent.replace_children([child, exchanged_child])
    return swapped_tree.refresh()


def _prefer_scored_move(
    left: _ScoredNniMove,
    right: _ScoredNniMove,
) -> bool:
    if left.score < right.score and not math.isclose(left.score, right.score):
        return True
    if right.score < left.score and not math.isclose(left.score, right.score):
        return False
    return left.newick < right.newick


def _node_sort_key(node: TreeNode) -> tuple[int, tuple[str, ...]]:
    descendants = tuple(descendant_taxa(node))
    return (len(descendants), descendants)


def _clade_id(node: TreeNode) -> str:
    return canonical_clade_id(frozenset(descendant_taxa(node)))


def _require_node_id(node: TreeNode) -> str:
    if node.node_id is None:
        raise AssertionError("NNI search requires refreshed node identities")
    return node.node_id
