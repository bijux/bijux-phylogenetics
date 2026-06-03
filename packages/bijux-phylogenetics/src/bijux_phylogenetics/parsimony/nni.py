from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.topology.rooted_nni import (
    RootedNniMoveCandidate,
    apply_rooted_nni_move,
    iter_rooted_nni_move_candidates,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .models import (
    FitchCharacterMatrix,
    ParsimonyCharacterWeights,
    ParsimonyNniSearchReport,
    ParsimonyNniTraceRow,
    SankoffCostMatrix,
)
from .topology_search import (
    resolve_topology_search_cost_matrix,
    resolve_topology_search_matrix,
    resolve_topology_search_method,
    resolve_topology_search_tree,
    resolve_topology_search_weights,
    score_topology_search_tree,
    topology_search_prefer_score,
    validate_topology_search_tree,
)


@dataclass(frozen=True, slots=True)
class _ScoredNniMove:
    candidate: RootedNniMoveCandidate
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
    character_weights: ParsimonyCharacterWeights
    | Mapping[str, float]
    | Path
    | None = None,
) -> ParsimonyNniSearchReport:
    """Search one rooted binary tree by accepting score-improving rooted NNI moves."""
    resolved_tree, resolved_tree_path = resolve_topology_search_tree(tree)
    validate_topology_search_tree(
        resolved_tree,
        workflow_name="parsimony NNI search",
    )
    resolved_matrix = resolve_topology_search_matrix(
        matrix,
        taxon_column=taxon_column,
    )
    resolved_method = resolve_topology_search_method(
        method,
        workflow_name="parsimony NNI search",
    )
    resolved_cost_matrix = resolve_topology_search_cost_matrix(
        method=resolved_method,
        cost_matrix=cost_matrix,
        allow_asymmetric_costs=allow_asymmetric_costs,
        workflow_name="parsimony NNI search",
    )
    resolved_weights = resolve_topology_search_weights(
        resolved_matrix,
        character_weights,
    )
    current_tree = resolved_tree.copy().refresh()
    current_score = score_topology_search_tree(
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
        for candidate in iter_rooted_nni_move_candidates(current_tree):
            neighbor_tree = apply_rooted_nni_move(current_tree, candidate)
            neighbor_score = score_topology_search_tree(
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
            if improving_move is None or topology_search_prefer_score(
                scored_move.score,
                scored_move.newick,
                improving_move.score,
                improving_move.newick,
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
