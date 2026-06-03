from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.trees import (
    analyze_branch_length_distribution,
    analyze_tree_set_branch_lengths,
)


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_analyze_branch_length_distribution_reports_summary_and_rows() -> None:
    report = analyze_branch_length_distribution(tree_fixture("example_tree.nwk"))

    assert report.tree_count == 1
    assert report.aggregate.branch_count == 6
    assert report.aggregate.defined_branch_count == 6
    assert report.aggregate.missing_branch_count == 0
    assert report.aggregate.zero_length_branch_count == 0
    assert report.aggregate.long_outlier_count == 0
    assert report.aggregate.minimum_branch_length == 0.1
    assert report.aggregate.maximum_branch_length == 0.2
    assert report.aggregate.mean_branch_length == 0.15
    assert report.aggregate.median_branch_length == 0.15
    assert report.aggregate.positive_branch_median == 0.15

    by_node = {row.node: row for row in report.rows}
    assert by_node["A"].branch_type == "terminal"
    assert by_node["A"].tip_taxon == "A"
    assert by_node["A"].descendant_tip_count == 1
    assert by_node["A"].root_depth == 0.3
    assert by_node["A"].outlier_class == "typical"
    assert by_node["A|B"].branch_type == "internal"
    assert by_node["A|B"].tip_taxon is None
    assert by_node["A|B"].descendant_tip_count == 2
    assert by_node["A|B"].root_depth == 0.2


def test_analyze_branch_length_distribution_flags_zero_negative_and_long_outliers() -> (
    None
):
    long_report = analyze_branch_length_distribution(
        tree_fixture("example_tree_long_branch.nwk")
    )
    zero_report = analyze_branch_length_distribution(
        tree_fixture("example_tree_zero_lengths.nwk")
    )
    negative_report = analyze_branch_length_distribution(
        tree_fixture("example_tree_negative_length.nwk")
    )

    assert long_report.aggregate.long_outlier_count == 1
    assert [row.node for row in long_report.rows if row.long_outlier] == ["A"]

    assert zero_report.aggregate.zero_length_branch_count == 3
    assert sorted(row.node for row in zero_report.rows if row.zero_length) == [
        "A",
        "A|B",
        "D",
    ]
    assert zero_report.aggregate.median_branch_length == 0.05

    assert negative_report.aggregate.negative_branch_count == 1
    assert [row.node for row in negative_report.rows if row.negative_length] == ["A"]


def test_analyze_tree_set_branch_lengths_aggregates_across_tree_rows() -> None:
    report = analyze_tree_set_branch_lengths(
        tree_fixture("example_tree_set_branch_lengths.nwk")
    )

    assert report.tree_count == 3
    assert report.aggregate.branch_count == 18
    assert report.aggregate.defined_branch_count == 18
    assert report.aggregate.zero_length_branch_count == 3
    assert report.aggregate.long_outlier_count == 1
    assert report.aggregate.minimum_branch_length == 0.0
    assert report.aggregate.maximum_branch_length == 1.0
    assert sorted({row.tree_index for row in report.rows}) == [1, 2, 3]
    assert [row.node for row in report.rows if row.long_outlier] == ["A"]
    assert [row.tree_index for row in report.rows if row.long_outlier] == [2]
