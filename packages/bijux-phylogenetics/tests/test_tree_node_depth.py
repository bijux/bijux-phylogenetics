from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.phylo.branch_lengths.node_depths import (
    compute_tree_node_depths,
    write_tree_node_depth_table,
)
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology import root_tree_on_outgroup
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_compute_tree_node_depths_matches_rooted_ultrametric_tree() -> None:
    report = compute_tree_node_depths(fixture("example_tree.nwk"))

    assert report.tip_labels == ["A", "B", "C", "D"]
    assert report.rooted is True
    assert report.tree_is_ultrametric is True
    assert report.node_count == 7
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.zero_branch_length_count == 0
    row_lookup = {row.node_id: row for row in report.rows}
    assert math.isclose(row_lookup[1].branch_length_depth, 0.3, abs_tol=1e-12)
    assert math.isclose(row_lookup[2].branch_length_depth, 0.3, abs_tol=1e-12)
    assert math.isclose(row_lookup[5].branch_length_depth, 0.0, abs_tol=1e-12)
    assert row_lookup[5].node_kind == "root"
    assert row_lookup[6].descendant_taxa == ["A", "B"]
    assert math.isclose(row_lookup[6].branch_length_depth, 0.2, abs_tol=1e-12)
    assert row_lookup[7].descendant_taxa == ["C", "D"]
    assert math.isclose(row_lookup[7].branch_length_depth, 0.1, abs_tol=1e-12)


def test_compute_tree_node_depths_handles_non_ultrametric_tree() -> None:
    report = compute_tree_node_depths(fixture("example_tree_ladderized.nwk"))

    assert report.rooted is True
    assert report.tree_is_ultrametric is False
    assert math.isclose(report.minimum_tip_depth, 0.1, abs_tol=1e-12)
    assert math.isclose(report.maximum_tip_depth, 0.3, abs_tol=1e-12)
    row_lookup = {row.node_id: row for row in report.rows}
    assert math.isclose(row_lookup[4].branch_length_depth, 0.1, abs_tol=1e-12)
    assert math.isclose(row_lookup[6].branch_length_depth, 0.1, abs_tol=1e-12)
    assert row_lookup[5].descendant_taxa == ["A", "B", "C", "D"]


def test_compute_tree_node_depths_handles_zero_branch_lengths() -> None:
    report = compute_tree_node_depths(fixture("example_tree_zero_lengths.nwk"))

    assert report.zero_branch_length_count > 0
    assert report.minimum_tip_depth == 0.0
    assert report.minimum_internal_depth == 0.0


def test_compute_tree_node_depths_rejects_missing_branch_lengths() -> None:
    with pytest.raises(InvalidBranchLengthError) as error:
        compute_tree_node_depths(fixture("example_tree_branch_lengths_missing.nwk"))

    assert error.value.code == "tree_node_depth_missing_branch_lengths"


def test_compute_tree_node_depths_after_outgroup_rooting(tmp_path: Path) -> None:
    rooted_tree, _report = root_tree_on_outgroup(
        fixture("example_tree.nwk"),
        outgroup_taxa=["D"],
    )
    rooted_path = tmp_path / "rooted-on-d.nwk"
    write_newick(rooted_path, rooted_tree)

    report = compute_tree_node_depths(rooted_path)

    assert report.rooted is True
    assert report.node_count == 7
    assert report.rows[4].node_kind == "root"
    assert report.rows[4].descendant_taxa == ["A", "B", "C", "D"]


def test_compute_tree_node_depths_after_pruning(tmp_path: Path) -> None:
    pruned_tree, _report = prune_tree_to_requested_taxa(
        fixture("example_tree_rooted_on_d.nwk"),
        ["A", "B", "C"],
    )
    pruned_path = tmp_path / "pruned-abc.nwk"
    write_newick(pruned_path, pruned_tree)

    report = compute_tree_node_depths(pruned_path)

    assert report.rooted is True
    assert report.tip_labels == ["A", "B", "C"]
    assert [row.node_id for row in report.rows] == [1, 2, 3, 4, 5]
    assert report.rows[3].node_kind == "root"
    assert report.rows[4].descendant_taxa == ["A", "B"]


def test_write_tree_node_depth_table_preserves_ape_node_order(tmp_path: Path) -> None:
    report = compute_tree_node_depths(fixture("example_tree.nwk"))
    output_path = tmp_path / "node-depths.tsv"

    write_tree_node_depth_table(output_path, report)

    rows = output_path.read_text(encoding="utf-8").splitlines()
    assert rows[0].startswith("node_id\tnode_kind\tnode_label\tdescendant_taxa")
    assert rows[1].startswith("1\ttip\tA\tA\t0.3")
    assert rows[5].startswith("5\troot\t\tA|B|C|D\t0")
