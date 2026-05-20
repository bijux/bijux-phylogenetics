from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.branch_lengths.branching_times import (
    compute_tree_branching_times,
    write_tree_branching_time_table,
)
from bijux_phylogenetics.runtime.errors import NonUltrametricTreeError

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


def test_compute_tree_branching_times_matches_small_ultrametric_tree() -> None:
    report = compute_tree_branching_times(fixture("example_tree.nwk"))

    assert report.rooted is True
    assert report.tree_is_ultrametric is True
    assert math.isclose(report.root_age, 0.3, abs_tol=1e-12)
    assert report.internal_node_count == 3
    row_lookup = {row.node_id: row for row in report.rows}
    assert math.isclose(row_lookup[5].branching_time, 0.3, abs_tol=1e-12)
    assert math.isclose(row_lookup[6].branching_time, 0.1, abs_tol=1e-12)
    assert math.isclose(row_lookup[7].branching_time, 0.2, abs_tol=1e-12)


def test_compute_tree_branching_times_preserves_internal_node_labels() -> None:
    report = compute_tree_branching_times(fixture("example_tree_named_clades.nwk"))

    row_lookup = {row.node_id: row for row in report.rows}
    assert row_lookup[5].node_label == "Root"
    assert row_lookup[6].node_label == "Mammals"
    assert row_lookup[7].node_label == "Birds"
    assert math.isclose(row_lookup[5].branching_time, 0.3, abs_tol=1e-12)
    assert math.isclose(row_lookup[6].branching_time, 0.1, abs_tol=1e-12)
    assert math.isclose(row_lookup[7].branching_time, 0.2, abs_tol=1e-12)


def test_compute_tree_branching_times_matches_medium_ultrametric_tree() -> None:
    report = compute_tree_branching_times(fixture("example_tree_eight_taxa.nwk"))

    assert math.isclose(report.root_age, 3.0, abs_tol=1e-12)
    assert report.internal_node_count == 7
    row_lookup = {row.node_id: row for row in report.rows}
    assert math.isclose(row_lookup[9].branching_time, 3.0, abs_tol=1e-12)
    assert math.isclose(row_lookup[10].branching_time, 2.0, abs_tol=1e-12)
    assert math.isclose(row_lookup[11].branching_time, 1.0, abs_tol=1e-12)
    assert math.isclose(row_lookup[15].branching_time, 1.0, abs_tol=1e-12)


def test_compute_tree_branching_times_handles_zero_internal_branch() -> None:
    report = compute_tree_branching_times(
        fixture("example_tree_ultrametric_zero_internal.nwk")
    )

    assert math.isclose(report.root_age, 0.3, abs_tol=1e-12)
    assert report.zero_branch_length_count == 1
    row_lookup = {row.node_id: row for row in report.rows}
    assert row_lookup[7].node_label == "Birds"
    assert math.isclose(row_lookup[7].branching_time, 0.3, abs_tol=1e-12)


def test_compute_tree_branching_times_rejects_non_ultrametric_tree() -> None:
    with pytest.raises(NonUltrametricTreeError) as error:
        compute_tree_branching_times(fixture("example_tree_ladderized.nwk"))

    assert error.value.code == "tree_branching_times_require_ultrametric_tree"
    assert error.value.details["max_tip_depth_deviation"] == pytest.approx(0.2)
    assert error.value.details["offending_taxa"] == ["A", "B", "D"]


def test_write_tree_branching_time_table_preserves_internal_node_order(
    tmp_path: Path,
) -> None:
    report = compute_tree_branching_times(fixture("example_tree_named_clades.nwk"))
    output_path = tmp_path / "branching-times.tsv"

    write_tree_branching_time_table(output_path, report)

    rows = output_path.read_text(encoding="utf-8").splitlines()
    assert rows[0].startswith("node_id\tnode_kind\tnode_label\tdescendant_taxa")
    assert rows[1].startswith("5\troot\tRoot\tA|B|C|D\t0\t0.3")
    assert rows[2].startswith("6\tinternal\tMammals\tA|B\t0.2\t0.1")
