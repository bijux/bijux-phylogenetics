from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.topology import prune_trees_to_shared_taxa


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_prune_trees_to_shared_taxa_preserves_pruning_audits_and_comparison() -> None:
    left, right, report = prune_trees_to_shared_taxa(
        fixture("example_tree.nwk"),
        fixture("example_tree_overlap.nwk"),
    )

    assert left.tip_names == ["A", "B", "C"]
    assert right.tip_names == ["A", "B", "C"]
    assert left.rooted is True
    assert right.rooted is True
    assert report.shared_taxa == ["A", "B", "C"]
    assert report.left_only_taxa == ["D"]
    assert report.right_only_taxa == ["E"]

    assert report.left_pruning.removed_taxa == ["D"]
    assert report.right_pruning.removed_taxa == ["E"]
    assert report.left_pruning.removed_taxa_with_reasons[0].reason == "not_requested"
    assert report.right_pruning.removed_taxa_with_reasons[0].reason == "not_requested"

    assert report.left_pruning.pruning_audit.root_to_tip_complete is True
    assert report.right_pruning.pruning_audit.root_to_tip_complete is True
    assert report.left_pruning.pruning_audit.pruned_total_branch_length == 0.7
    assert report.right_pruning.pruning_audit.pruned_total_branch_length == 0.7
    assert report.left_pruning.pruning_audit.branch_length_delta == -0.2
    assert report.right_pruning.pruning_audit.branch_length_delta == -0.2

    assert report.post_pruning_comparison.shared_taxa == ["A", "B", "C"]
    assert report.post_pruning_comparison.left_only_taxa == []
    assert report.post_pruning_comparison.right_only_taxa == []
    assert report.post_pruning_comparison.taxon_overlap_policy == "require-identical"
    assert report.post_pruning_comparison.robinson_foulds_distance == 0
    assert report.post_pruning_comparison.topology_equal is True
    assert (
        report.post_pruning_comparison.same_topology_different_branch_lengths is False
    )
