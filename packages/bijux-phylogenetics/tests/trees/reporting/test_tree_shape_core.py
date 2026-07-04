from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    summarize_tree_set_shapes,
    summarize_tree_shape,
)


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_summarize_tree_shape_reports_balanced_tree_metrics() -> None:
    report = summarize_tree_shape(tree_fixture("example_tree.nwk"))

    assert report.tree_count == 1
    assert report.source_format == "newick"
    row = report.rows[0]
    assert row.tree_index is None
    assert row.tip_count == 4
    assert row.internal_node_count == 3
    assert row.is_binary is True
    assert row.cherry_count == 2
    assert row.sackin_imbalance_index == 8
    assert row.colless_imbalance_index == 0.0
    assert row.normalized_colless_imbalance == 0.0
    assert row.tree_height_edges == 2
    assert row.tree_height_branch_length == 0.3
    assert row.mean_tip_depth_edges == 2.0
    assert row.mean_root_to_tip_branch_length == 0.3
    assert row.imbalance_summary == "balanced"
    assert row.ladderized is False
    assert row.star_like is False
    assert row.comb_like is False
    assert row.unusually_imbalanced is False

    assert report.aggregate.balanced_tree_count == 1
    assert report.aggregate.mean_sackin_imbalance_index == 8.0
    assert report.aggregate.mean_tree_height_branch_length == 0.3


def test_summarize_tree_shape_reports_ladderized_tree_metrics() -> None:
    report = summarize_tree_shape(tree_fixture("example_tree_ladderized.nwk"))

    row = report.rows[0]
    assert row.cherry_count == 1
    assert row.sackin_imbalance_index == 9
    assert row.colless_imbalance_index == 3.0
    assert row.normalized_colless_imbalance == 1.0
    assert row.tree_height_edges == 3
    assert row.tree_height_branch_length == 0.3
    assert row.mean_tip_depth_edges == 2.25
    assert row.mean_root_to_tip_branch_length == 0.225
    assert row.imbalance_summary == "ladderized"
    assert row.ladderized is True
    assert row.comb_like is True
    assert row.unusually_imbalanced is True


def test_summarize_tree_shape_keeps_colless_unavailable_for_polytomy() -> None:
    report = summarize_tree_shape(tree_fixture("example_tree_polytomy.nwk"))

    row = report.rows[0]
    assert row.is_binary is False
    assert row.colless_imbalance_index is None
    assert row.normalized_colless_imbalance is None
    assert row.unusually_imbalanced is None

    assert report.aggregate.colless_defined_tree_count == 0
    assert report.aggregate.mean_colless_imbalance_index is None


def test_summarize_tree_set_shapes_reports_one_row_per_tree_and_aggregate() -> None:
    report = summarize_tree_set_shapes(tree_fixture("example_tree_set_left.nwk"))

    assert report.tree_count == 3
    assert [row.tree_index for row in report.rows] == [1, 2, 3]
    assert all(row.imbalance_summary == "balanced" for row in report.rows)
    assert all(row.tree_height_edges == 2 for row in report.rows)
    assert report.aggregate.balanced_tree_count == 3
    assert report.aggregate.ladderized_tree_count == 0
    assert report.aggregate.mean_cherry_count == 2.0
    assert report.aggregate.mean_sackin_imbalance_index == 8.0
    assert report.aggregate.maximum_tree_height_edges == 2
