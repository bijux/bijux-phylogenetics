from __future__ import annotations

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.topology import (
    RootedSprEnumerationBudget,
    enumerate_rooted_spr_neighbors,
    iter_rooted_spr_move_candidates,
)


def test_rooted_spr_prune_budget_reports_skipped_pruned_clades() -> None:
    budget = RootedSprEnumerationBudget(max_pruned_clade_count=1)

    report = enumerate_rooted_spr_neighbors(loads_newick("(((A,C),B),D);"), budget=budget)

    assert report.max_pruned_clade_count == 1
    assert report.max_regraft_target_count_per_pruned_clade is None
    assert report.skipped_pruned_clade_count == 5
    assert report.skipped_regraft_target_count == 0
    assert report.generated_move_candidate_count == 6
    assert report.identity_move_candidate_count == 1
    assert report.self_regraft_candidate_count == 0
    assert report.generated_neighbor_count == 4
    assert report.unique_neighbor_topology_count == 4
    assert len(report.duplicate_move_neighbor_topologies) == 1


def test_rooted_spr_regraft_budget_reports_skipped_regraft_targets() -> None:
    budget = RootedSprEnumerationBudget(max_regraft_target_count_per_pruned_clade=3)

    report = enumerate_rooted_spr_neighbors(loads_newick("(((A,C),B),D);"), budget=budget)

    assert report.max_pruned_clade_count is None
    assert report.max_regraft_target_count_per_pruned_clade == 3
    assert report.skipped_pruned_clade_count == 0
    assert report.skipped_regraft_target_count == 13
    assert report.generated_move_candidate_count == 17
    assert report.identity_move_candidate_count == 6
    assert report.self_regraft_candidate_count == 0
    assert report.generated_neighbor_count == 7
    assert report.unique_neighbor_topology_count == 7
    assert len(report.duplicate_move_neighbor_topologies) == 4


def test_rooted_spr_combined_budget_limits_candidate_iterator_deterministically() -> None:
    budget = RootedSprEnumerationBudget(
        max_pruned_clade_count=1,
        max_regraft_target_count_per_pruned_clade=3,
    )

    candidates = list(
        iter_rooted_spr_move_candidates(loads_newick("(((A,C),B),D);"), budget=budget)
    )
    report = enumerate_rooted_spr_neighbors(loads_newick("(((A,C),B),D);"), budget=budget)

    assert len(candidates) == 3
    assert report.skipped_pruned_clade_count == 5
    assert report.skipped_regraft_target_count == 3
    assert report.generated_move_candidate_count == 3
    assert report.identity_move_candidate_count == 1
    assert report.self_regraft_candidate_count == 0
    assert report.generated_neighbor_count == 2
    assert report.unique_neighbor_topology_count == 2
    assert report.duplicate_move_neighbor_topologies == []

