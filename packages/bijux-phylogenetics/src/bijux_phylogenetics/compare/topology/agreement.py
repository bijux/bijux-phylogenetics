from __future__ import annotations

from collections.abc import Iterable
from itertools import combinations
import math
from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import _load_tree
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
    RobinsonFouldsMode,
)

AGREEMENT_SUBTREE_SEARCH_STRATEGY = "exact-descending-retained-subsets"


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


def _build_agreement_candidate_row(
    left: PhyloTree,
    right: PhyloTree,
    *,
    shared_taxa: list[str],
    retained_taxa: list[str],
    candidate_index: int,
    rf_mode: RobinsonFouldsMode,
) -> tuple[AgreementSubtreeCandidateRow, bool]:
    pruned_left = prune_tree_object_to_requested_taxa(left, retained_taxa)
    pruned_right = prune_tree_object_to_requested_taxa(right, retained_taxa)
    comparison = _compare_tree_objects(
        pruned_left,
        pruned_right,
        rf_mode=rf_mode,
        taxon_overlap_policy="require-identical",
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
