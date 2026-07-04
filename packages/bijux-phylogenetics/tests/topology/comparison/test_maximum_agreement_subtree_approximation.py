from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.compare as compare_api
from bijux_phylogenetics.compare import (
    MaximumAgreementSubtreeApproximationReport,
    MaximumAgreementSubtreeSearchRow,
    approximate_maximum_agreement_subtree,
    prune_trees_to_agreement_subtree,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_package_compare_gateway_exports_maximum_agreement_subtree_surface() -> None:
    assert (
        compare_api.MaximumAgreementSubtreeApproximationReport
        is MaximumAgreementSubtreeApproximationReport
    )
    assert (
        compare_api.MaximumAgreementSubtreeSearchRow is MaximumAgreementSubtreeSearchRow
    )
    assert (
        compare_api.approximate_maximum_agreement_subtree
        is approximate_maximum_agreement_subtree
    )


def test_approximate_maximum_agreement_subtree_matches_exact_small_fixture() -> None:
    exact_left, exact_right, exact_report = prune_trees_to_agreement_subtree(
        fixture("agreement_subtree_left.nwk"),
        fixture("agreement_subtree_right.nwk"),
    )
    left, right, report = approximate_maximum_agreement_subtree(
        fixture("agreement_subtree_left.nwk"),
        fixture("agreement_subtree_right.nwk"),
        max_evaluated_candidate_count=13,
    )

    assert left.tip_names == ["A", "B", "D", "E"]
    assert right.tip_names == ["A", "B", "D", "E"]
    assert left.to_newick() == exact_left.to_newick()
    assert right.to_newick() == exact_right.to_newick()
    assert report.shared_taxa == ["A", "B", "C", "D", "E"]
    assert report.retained_taxa == exact_report.retained_taxa == ["A", "B", "D", "E"]
    assert report.approximation_removed_taxa == ["C"]
    assert report.search_strategy == "greedy-single-taxon-removal"
    assert (
        report.selection_objective
        == "minimize-robinson-foulds-then-normalized-distance"
    )
    assert report.approximation_status == "heuristic-solution-not-guaranteed-optimal"
    assert report.possible_retained_subset_count == 26
    assert report.max_evaluated_candidate_count == 13
    assert report.evaluated_candidate_count == 6
    assert [
        row.retained_taxa for row in report.search_rows if row.selected_for_next_step
    ] == [
        ["A", "B", "C", "D", "E"],
        ["A", "B", "D", "E"],
    ]
    assert report.post_pruning_comparison.robinson_foulds_distance == 0
    assert report.post_pruning_comparison.topology_equal is True


def test_approximate_maximum_agreement_subtree_requires_positive_budget() -> None:
    with pytest.raises(ValueError, match="positive explicit candidate budget"):
        approximate_maximum_agreement_subtree(
            fixture("agreement_subtree_left.nwk"),
            fixture("agreement_subtree_right.nwk"),
            max_evaluated_candidate_count=0,
        )


def test_approximate_maximum_agreement_subtree_rejects_too_small_frontier_budget() -> (
    None
):
    with pytest.raises(ValueError, match="need at least 6 candidates"):
        approximate_maximum_agreement_subtree(
            fixture("agreement_subtree_left.nwk"),
            fixture("agreement_subtree_right.nwk"),
            max_evaluated_candidate_count=5,
        )


def test_approximate_maximum_agreement_subtree_reports_large_fixture_status() -> None:
    _left, _right, report = approximate_maximum_agreement_subtree(
        fixture("maximum_agreement_large_left.nwk"),
        fixture("maximum_agreement_large_right.nwk"),
        max_evaluated_candidate_count=32,
    )

    assert report.shared_taxa == ["A", "B", "C", "D", "E", "F", "G", "H"]
    assert report.approximation_status == "heuristic-solution-not-guaranteed-optimal"
    assert report.max_evaluated_candidate_count == 32
    assert report.evaluated_candidate_count <= 32
    assert len(report.retained_taxa) < len(report.shared_taxa)
    assert report.post_pruning_comparison.robinson_foulds_distance == 0
    assert report.post_pruning_comparison.topology_equal is True
