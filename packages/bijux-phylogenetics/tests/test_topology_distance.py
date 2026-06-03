from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.compare.topology import (
    compare_topology_distance_trees,
    write_topology_distance_split_table,
)
from bijux_phylogenetics.datasets.shared_fixtures import get_shared_tree_set_fixture
from bijux_phylogenetics.io.newick import load_newick_tree_set, write_newick


def _fixture_pair(tmp_path: Path, fixture_id: str) -> tuple[Path, Path]:
    fixture = get_shared_tree_set_fixture(fixture_id)
    trees = load_newick_tree_set(fixture.path)
    assert len(trees) == 2
    left_path = tmp_path / f"{fixture_id}-left.nwk"
    right_path = tmp_path / f"{fixture_id}-right.nwk"
    write_newick(left_path, trees[0])
    write_newick(right_path, trees[1])
    return left_path, right_path


def test_compare_topology_distance_reports_zero_distance_for_identical_rooted_pair(
    tmp_path: Path,
) -> None:
    left_path, right_path = _fixture_pair(
        tmp_path, "topology_distance_identical_rooted_pair"
    )
    left_tree, right_tree = load_newick_tree_set(
        get_shared_tree_set_fixture("topology_distance_identical_rooted_pair").path
    )

    report = compare_topology_distance_trees(
        left_tree,
        right_tree,
        left_path=left_path,
        right_path=right_path,
    )

    assert report.topology_equal is True
    assert report.robinson_foulds_distance == 0
    assert report.left_split_count == 2
    assert report.shared_split_count == 2
    assert report.left_only_split_count == 0
    assert report.right_only_split_count == 0
    assert [row.comparison_status for row in report.split_rows] == ["shared", "shared"]


def test_compare_topology_distance_ignores_child_order_for_rooted_pair(
    tmp_path: Path,
) -> None:
    left_path, right_path = _fixture_pair(
        tmp_path, "topology_distance_rooted_child_order_pair"
    )
    left_tree, right_tree = load_newick_tree_set(
        get_shared_tree_set_fixture("topology_distance_rooted_child_order_pair").path
    )

    report = compare_topology_distance_trees(
        left_tree,
        right_tree,
        left_path=left_path,
        right_path=right_path,
    )

    assert report.topology_equal is True
    assert report.robinson_foulds_distance == 0
    assert report.left_split_count == report.right_split_count == 2


def test_compare_topology_distance_reports_rooted_conflict_rows(tmp_path: Path) -> None:
    left_path, right_path = _fixture_pair(
        tmp_path, "topology_distance_rooted_conflict_pair"
    )
    left_tree, right_tree = load_newick_tree_set(
        get_shared_tree_set_fixture("topology_distance_rooted_conflict_pair").path
    )

    report = compare_topology_distance_trees(
        left_tree,
        right_tree,
        left_path=left_path,
        right_path=right_path,
    )

    assert report.topology_equal is False
    assert report.robinson_foulds_distance == 2
    assert report.normalized_robinson_foulds == pytest.approx(0.5)
    assert report.shared_split_count == 1
    assert report.left_only_split_count == 1
    assert report.right_only_split_count == 1
    assert [(row.split_id, row.comparison_status) for row in report.split_rows] == [
        ("A|B", "shared"),
        ("C|D", "left_only"),
        ("A|B|C", "right_only"),
    ]


def test_compare_topology_distance_reports_polytomy_presence(tmp_path: Path) -> None:
    left_path, right_path = _fixture_pair(
        tmp_path, "topology_distance_rooted_polytomy_pair"
    )
    left_tree, right_tree = load_newick_tree_set(
        get_shared_tree_set_fixture("topology_distance_rooted_polytomy_pair").path
    )

    report = compare_topology_distance_trees(
        left_tree,
        right_tree,
        left_path=left_path,
        right_path=right_path,
    )

    assert report.polytomy_present_left is False
    assert report.polytomy_present_right is True
    assert report.robinson_foulds_distance == 3
    assert report.left_only_split_count == 2
    assert report.right_only_split_count == 1


def test_compare_topology_distance_supports_unrooted_pairs(tmp_path: Path) -> None:
    left_path, right_path = _fixture_pair(
        tmp_path, "topology_distance_unrooted_conflict_pair"
    )
    left_tree, right_tree = load_newick_tree_set(
        get_shared_tree_set_fixture("topology_distance_unrooted_conflict_pair").path
    )

    report = compare_topology_distance_trees(
        left_tree,
        right_tree,
        left_path=left_path,
        right_path=right_path,
        rf_mode="unrooted",
    )

    assert report.rf_mode == "unrooted"
    assert report.robinson_foulds_distance == 2
    assert report.normalized_robinson_foulds == pytest.approx(1.0)
    assert report.left_split_count == report.right_split_count == 1
    assert all(row.split_kind == "split" for row in report.split_rows)


def test_compare_topology_distance_scales_to_large_rooted_pair(tmp_path: Path) -> None:
    left_path, right_path = _fixture_pair(
        tmp_path, "topology_distance_large_rooted_pair"
    )
    left_tree, right_tree = load_newick_tree_set(
        get_shared_tree_set_fixture("topology_distance_large_rooted_pair").path
    )

    report = compare_topology_distance_trees(
        left_tree,
        right_tree,
        left_path=left_path,
        right_path=right_path,
    )

    assert len(report.shared_taxa) == 128
    assert report.left_split_count == report.right_split_count == 126
    assert report.robinson_foulds_distance == 24
    assert report.normalized_robinson_foulds == pytest.approx(24 / 252)


def test_write_topology_distance_split_table_writes_governed_rows(
    tmp_path: Path,
) -> None:
    left_path, right_path = _fixture_pair(
        tmp_path, "topology_distance_rooted_conflict_pair"
    )
    output_path = tmp_path / "topology-distance.tsv"

    write_topology_distance_split_table(output_path, left_path, right_path)

    assert output_path.exists()
    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "split_id\tsplit_kind\tcomparison_status\ttaxon_count\tdescendant_taxa\tleft_present\tright_present",
        "A|B\tclade\tshared\t2\tA|B\ttrue\ttrue",
        "C|D\tclade\tleft_only\t2\tC|D\ttrue\tfalse",
        "A|B|C\tclade\tright_only\t3\tA|B|C\tfalse\ttrue",
    ]
