from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.compare.topology import compare_topology_distance_trees
from bijux_phylogenetics.datasets.shared_fixtures import get_shared_tree_set_fixture
from bijux_phylogenetics.io.newick import load_newick_tree_set, write_newick
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

pytestmark = pytest.mark.smoke


def test_native_tree_model_newick_helpers_roundtrip_deterministically() -> None:
    tree = PhyloTree.from_newick("((B:0.1,C:0.1)Mammals:0.2,A:0.1)Root;")

    assert tree.to_newick() == "(A:0.1,(B:0.1,C:0.1)Mammals:0.2)Root;"
    assert tree.root.name == "Root"
    assert [node.name for node in tree.iter_internal_nodes()] == ["Root", "Mammals"]


def test_topology_distance_reports_rooted_conflict_rows(tmp_path: Path) -> None:
    fixture = get_shared_tree_set_fixture("topology_distance_rooted_conflict_pair")
    left_tree, right_tree = load_newick_tree_set(fixture.path)
    left_path = tmp_path / "left.nwk"
    right_path = tmp_path / "right.nwk"
    write_newick(left_path, left_tree)
    write_newick(right_path, right_tree)

    report = compare_topology_distance_trees(
        left_tree,
        right_tree,
        left_path=left_path,
        right_path=right_path,
    )

    assert report.topology_equal is False
    assert report.robinson_foulds_distance == 2
    assert report.left_only_split_count == 1
    assert report.right_only_split_count == 1
