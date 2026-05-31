from __future__ import annotations

from collections.abc import Iterable
from itertools import combinations
import math
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation.structure import _load_tree
from bijux_phylogenetics.phylo.pruning import (
    prune_tree_object_to_requested_taxa,
    prune_tree_to_requested_taxa,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .comparison import (
    _build_tree_comparison_report,
    _compare_tree_objects,
    _resolve_shared_taxa,
)
from .models import (
    AgreementSubtreeCandidateRow,
    AgreementSubtreePruningReport,
    MaximumAgreementSubtreeApproximationReport,
    MaximumAgreementSubtreeSearchRow,
    RobinsonFouldsMode,
)

AGREEMENT_SUBTREE_SEARCH_STRATEGY = "exact-descending-retained-subsets"
MAXIMUM_AGREEMENT_SUBTREE_SEARCH_STRATEGY = "greedy-single-taxon-removal"
MAXIMUM_AGREEMENT_SUBTREE_SELECTION_OBJECTIVE = (
    "minimize-robinson-foulds-then-normalized-distance"
)
MAXIMUM_AGREEMENT_SUBTREE_APPROXIMATION_STATUS = (
    "heuristic-solution-not-guaranteed-optimal"
)


def _candidate_subset_count(shared_taxon_count: int) -> int:
    return sum(
        math.comb(shared_taxon_count, retained)
        for retained in range(2, shared_taxon_count + 1)
    )


def _iter_retained_taxon_subsets(shared_taxa: list[str]) -> Iterable[list[str]]:
    """Yield largest retained taxon subsets first with lexicographic tie-breaking."""
    ordered_taxa = sorted(shared_taxa)
    for retained_taxon_count in range(len(ordered_taxa), 1, -1):
        for retained_taxa in combinations(ordered_taxa, retained_taxon_count):
            yield list(retained_taxa)


def _evaluate_agreement_subset(
    left: PhyloTree,
    right: PhyloTree,
    *,
    retained_taxa: list[str],
    rf_mode: RobinsonFouldsMode,
):
    pruned_left = prune_tree_object_to_requested_taxa(left, retained_taxa)
    pruned_right = prune_tree_object_to_requested_taxa(right, retained_taxa)
    return _compare_tree_objects(
        pruned_left,
        pruned_right,
        rf_mode=rf_mode,
        taxon_overlap_policy="require-identical",
    )


def _build_agreement_candidate_row(
    left: PhyloTree,
    right: PhyloTree,
    *,
    shared_taxa: list[str],
    retained_taxa: list[str],
    candidate_index: int,
    rf_mode: RobinsonFouldsMode,
) -> tuple[AgreementSubtreeCandidateRow, bool]:
    comparison = _evaluate_agreement_subset(
        left,
        right,
        retained_taxa=retained_taxa,
        rf_mode=rf_mode,
    )
    topology_equal = comparison.robinson_foulds_distance == 0
    removed_taxa = sorted(set(shared_taxa) - set(retained_taxa))
    return (
        AgreementSubtreeCandidateRow(
            candidate_index=candidate_index,
            retained_taxon_count=len(retained_taxa),
            retained_taxa=list(retained_taxa),
            removed_taxa=removed_taxa,
            robinson_foulds_distance=comparison.robinson_foulds_distance,
            normalized_robinson_foulds=comparison.normalized_robinson_foulds,
            topology_equal=topology_equal,
        ),
        topology_equal,
    )


def prune_trees_to_agreement_subtree(
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
) -> tuple[PhyloTree, PhyloTree, AgreementSubtreePruningReport]:
    """Prune two conflicting trees to the largest exact shared agreement subtree."""
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    shared_taxa_set, left_only_taxa, right_only_taxa = _resolve_shared_taxa(
        set(left.tip_names),
        set(right.tip_names),
        taxon_overlap_policy="prune-to-shared",
    )
    shared_taxa = sorted(shared_taxa_set)
    candidate_rows: list[AgreementSubtreeCandidateRow] = []
    retained_taxa: list[str] | None = None
    for candidate_index, retained_subset in enumerate(
        _iter_retained_taxon_subsets(shared_taxa),
        start=1,
    ):
        candidate_row, topology_equal = _build_agreement_candidate_row(
            left,
            right,
            shared_taxa=shared_taxa,
            retained_taxa=retained_subset,
            candidate_index=candidate_index,
            rf_mode=rf_mode,
        )
        candidate_rows.append(candidate_row)
        if topology_equal:
            retained_taxa = retained_subset
            break
    if retained_taxa is None:
        raise AssertionError(
            "agreement subtree pruning must find a matching retained taxon subset"
        )

    pruned_left, left_pruning = prune_tree_to_requested_taxa(
        left_path,
        retained_taxa,
    )
    pruned_right, right_pruning = prune_tree_to_requested_taxa(
        right_path,
        retained_taxa,
    )
    post_pruning_comparison = _build_tree_comparison_report(
        left_path,
        right_path,
        pruned_left,
        pruned_right,
        rf_mode=rf_mode,
        taxon_overlap_policy="require-identical",
    )
    if post_pruning_comparison.robinson_foulds_distance != 0:
        raise AssertionError(
            "agreement subtree pruning must return matching pruned topologies"
        )
    report = AgreementSubtreePruningReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=shared_taxa,
        left_only_taxa=left_only_taxa,
        right_only_taxa=right_only_taxa,
        rf_mode=rf_mode,
        search_strategy=AGREEMENT_SUBTREE_SEARCH_STRATEGY,
        possible_retained_subset_count=_candidate_subset_count(len(shared_taxa)),
        evaluated_candidate_count=len(candidate_rows),
        retained_taxa=retained_taxa,
        agreement_removed_taxa=sorted(set(shared_taxa) - set(retained_taxa)),
        left_pruning=left_pruning,
        right_pruning=right_pruning,
        post_pruning_comparison=post_pruning_comparison,
        candidate_rows=candidate_rows,
    )
    return pruned_left, pruned_right, report


def _build_maximum_agreement_search_row(
    *,
    evaluation_index: int,
    step_index: int,
    shared_taxa: list[str],
    retained_taxa: list[str],
    robinson_foulds_distance: int,
    normalized_robinson_foulds: float,
    topology_equal: bool,
    selected_for_next_step: bool,
) -> MaximumAgreementSubtreeSearchRow:
    return MaximumAgreementSubtreeSearchRow(
        evaluation_index=evaluation_index,
        step_index=step_index,
        retained_taxon_count=len(retained_taxa),
        retained_taxa=list(retained_taxa),
        removed_taxa=sorted(set(shared_taxa) - set(retained_taxa)),
        robinson_foulds_distance=robinson_foulds_distance,
        normalized_robinson_foulds=normalized_robinson_foulds,
        topology_equal=topology_equal,
        selected_for_next_step=selected_for_next_step,
    )


def approximate_maximum_agreement_subtree(
    left_path: Path,
    right_path: Path,
    *,
    rf_mode: RobinsonFouldsMode = "rooted",
    max_evaluated_candidate_count: int,
) -> tuple[PhyloTree, PhyloTree, MaximumAgreementSubtreeApproximationReport]:
    """Approximate one large agreement subtree under an explicit candidate budget."""
    if max_evaluated_candidate_count <= 0:
        raise ValueError(
            "maximum agreement subtree approximation requires a positive explicit candidate budget"
        )
    left = _load_tree(left_path)
    right = _load_tree(right_path)
    shared_taxa_set, left_only_taxa, right_only_taxa = _resolve_shared_taxa(
        set(left.tip_names),
        set(right.tip_names),
        taxon_overlap_policy="prune-to-shared",
    )
    shared_taxa = sorted(shared_taxa_set)
    evaluated_candidate_count = 0
    current_retained_taxa = list(shared_taxa)
    current_comparison = _evaluate_agreement_subset(
        left,
        right,
        retained_taxa=current_retained_taxa,
        rf_mode=rf_mode,
    )
    evaluated_candidate_count += 1
    search_rows = [
        _build_maximum_agreement_search_row(
            evaluation_index=evaluated_candidate_count,
            step_index=0,
            shared_taxa=shared_taxa,
            retained_taxa=current_retained_taxa,
            robinson_foulds_distance=current_comparison.robinson_foulds_distance,
            normalized_robinson_foulds=current_comparison.normalized_robinson_foulds,
            topology_equal=current_comparison.robinson_foulds_distance == 0,
            selected_for_next_step=True,
        )
    ]
    step_index = 0
    while current_comparison.robinson_foulds_distance != 0:
        if len(current_retained_taxa) <= 2:
            raise AssertionError(
                "maximum agreement subtree approximation must reach agreement by two taxa"
            )
        next_candidate_count = len(current_retained_taxa)
        required_budget = evaluated_candidate_count + next_candidate_count
        if required_budget > max_evaluated_candidate_count:
            raise ValueError(
                "maximum agreement subtree approximation budget is too small "
                f"to evaluate the next removal frontier; need at least {required_budget} candidates"
            )
        step_index += 1
        step_rows: list[MaximumAgreementSubtreeSearchRow] = []
        best_rank: tuple[int, float, tuple[str, ...]] | None = None
        best_retained_taxa: list[str] | None = None
        best_comparison = None
        best_row_index = -1
        for taxon in current_retained_taxa:
            candidate_retained_taxa = [
                candidate_taxon
                for candidate_taxon in current_retained_taxa
                if candidate_taxon != taxon
            ]
            candidate_comparison = _evaluate_agreement_subset(
                left,
                right,
                retained_taxa=candidate_retained_taxa,
                rf_mode=rf_mode,
            )
            evaluated_candidate_count += 1
            step_rows.append(
                _build_maximum_agreement_search_row(
                    evaluation_index=evaluated_candidate_count,
                    step_index=step_index,
                    shared_taxa=shared_taxa,
                    retained_taxa=candidate_retained_taxa,
                    robinson_foulds_distance=(
                        candidate_comparison.robinson_foulds_distance
                    ),
                    normalized_robinson_foulds=(
                        candidate_comparison.normalized_robinson_foulds
                    ),
                    topology_equal=(candidate_comparison.robinson_foulds_distance == 0),
                    selected_for_next_step=False,
                )
            )
            candidate_rank = (
                candidate_comparison.robinson_foulds_distance,
                candidate_comparison.normalized_robinson_foulds,
                tuple(candidate_retained_taxa),
            )
            if best_rank is None or candidate_rank < best_rank:
                best_rank = candidate_rank
                best_retained_taxa = candidate_retained_taxa
                best_comparison = candidate_comparison
                best_row_index = len(step_rows) - 1
        if best_retained_taxa is None or best_comparison is None or best_row_index < 0:
            raise AssertionError(
                "maximum agreement subtree approximation must evaluate at least one child candidate"
            )
        step_rows[best_row_index].selected_for_next_step = True
        search_rows.extend(step_rows)
        current_retained_taxa = best_retained_taxa
        current_comparison = best_comparison

    pruned_left, left_pruning = prune_tree_to_requested_taxa(
        left_path,
        current_retained_taxa,
    )
    pruned_right, right_pruning = prune_tree_to_requested_taxa(
        right_path,
        current_retained_taxa,
    )
    post_pruning_comparison = _build_tree_comparison_report(
        left_path,
        right_path,
        pruned_left,
        pruned_right,
        rf_mode=rf_mode,
        taxon_overlap_policy="require-identical",
    )
    if post_pruning_comparison.robinson_foulds_distance != 0:
        raise AssertionError(
            "maximum agreement subtree approximation must return one agreement subtree"
        )
    report = MaximumAgreementSubtreeApproximationReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=shared_taxa,
        left_only_taxa=left_only_taxa,
        right_only_taxa=right_only_taxa,
        rf_mode=rf_mode,
        search_strategy=MAXIMUM_AGREEMENT_SUBTREE_SEARCH_STRATEGY,
        selection_objective=MAXIMUM_AGREEMENT_SUBTREE_SELECTION_OBJECTIVE,
        approximation_status=MAXIMUM_AGREEMENT_SUBTREE_APPROXIMATION_STATUS,
        possible_retained_subset_count=_candidate_subset_count(len(shared_taxa)),
        max_evaluated_candidate_count=max_evaluated_candidate_count,
        evaluated_candidate_count=evaluated_candidate_count,
        retained_taxa=current_retained_taxa,
        approximation_removed_taxa=sorted(
            set(shared_taxa) - set(current_retained_taxa)
        ),
        left_pruning=left_pruning,
        right_pruning=right_pruning,
        post_pruning_comparison=post_pruning_comparison,
        search_rows=search_rows,
    )
    return pruned_left, pruned_right, report
