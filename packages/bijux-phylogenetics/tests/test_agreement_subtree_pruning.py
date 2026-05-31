from __future__ import annotations

from pathlib import Path

import bijux_phylogenetics.compare as compare_api
from bijux_phylogenetics.compare import (
    AgreementSubtreeCandidateRow,
    AgreementSubtreePruningReport,
    prune_trees_to_agreement_subtree,
)
from bijux_phylogenetics.compare.topology import prune_trees_to_shared_taxa


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_package_compare_gateway_exports_agreement_subtree_surface() -> None:
    assert compare_api.AgreementSubtreeCandidateRow is AgreementSubtreeCandidateRow
    assert compare_api.AgreementSubtreePruningReport is AgreementSubtreePruningReport
    assert (
        compare_api.prune_trees_to_agreement_subtree is prune_trees_to_agreement_subtree
    )


def test_prune_trees_to_agreement_subtree_resolves_conflicting_shared_taxa() -> None:
    _shared_left, _shared_right, shared_report = prune_trees_to_shared_taxa(
        fixture("agreement_subtree_left.nwk"),
        fixture("agreement_subtree_right.nwk"),
    )
    left, right, report = prune_trees_to_agreement_subtree(
        fixture("agreement_subtree_left.nwk"),
        fixture("agreement_subtree_right.nwk"),
    )

    assert shared_report.shared_taxa == ["A", "B", "C", "D", "E"]
    assert shared_report.post_pruning_comparison.robinson_foulds_distance == 2
    assert left.tip_names == ["A", "B", "D", "E"]
    assert right.tip_names == ["A", "B", "D", "E"]
    assert left.to_newick() == "((A:0.1,B:0.1):0.2,(D:0.1,E:0.1):0.2);"
    assert right.to_newick() == "((A:0.1,B:0.1):0.2,(D:0.1,E:0.1):0.2);"
    assert report.shared_taxa == ["A", "B", "C", "D", "E"]
    assert report.retained_taxa == ["A", "B", "D", "E"]
    assert report.agreement_removed_taxa == ["C"]
    assert report.rf_mode == "rooted"
    assert report.search_strategy == "exact-descending-retained-subsets"
    assert report.possible_retained_subset_count == 26
    assert report.evaluated_candidate_count == 4
    assert report.left_pruning.removed_taxa == ["C"]
    assert report.right_pruning.removed_taxa == ["C"]
    assert report.post_pruning_comparison.robinson_foulds_distance == 0
    assert report.post_pruning_comparison.topology_equal is True
    assert [row.retained_taxa for row in report.candidate_rows] == [
        ["A", "B", "C", "D", "E"],
        ["A", "B", "C", "D"],
        ["A", "B", "C", "E"],
        ["A", "B", "D", "E"],
    ]
    assert [row.robinson_foulds_distance for row in report.candidate_rows] == [
        2,
        2,
        2,
        0,
    ]
    assert report.candidate_rows[-1].topology_equal is True
