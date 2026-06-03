from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.trees import extract_tree_clades, extract_tree_set_clades


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def metadata_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def test_extract_tree_clades_reports_support_depth_age_and_metadata() -> None:
    report = extract_tree_clades(
        tree_fixture("example_tree_support_conflict_left.nwk"),
        metadata_path=metadata_fixture("example_metadata.tsv"),
        metadata_columns=["species", "location"],
    )

    assert report.tree_count == 1
    assert report.metadata_columns == ["species", "location"]
    assert len(report.rows) == 7
    by_clade = {row.clade_id: row for row in report.rows}

    assert by_clade["A|B|C|D"].node_kind == "root"
    assert by_clade["A|B|C|D"].root_depth == 0.0
    assert by_clade["A|B|C|D"].descendant_tip_depth_max == 0.3
    assert by_clade["A|B|C|D"].node_age == 0.3

    assert by_clade["A|B"].node_kind == "internal"
    assert by_clade["A|B"].support == 95.0
    assert by_clade["A|B"].support_fraction == 0.95
    assert by_clade["A|B"].branch_length == 0.2
    assert by_clade["A|B"].root_depth == 0.2
    assert by_clade["A|B"].descendant_tip_depth_min == 0.1
    assert by_clade["A|B"].descendant_tip_depth_max == 0.1
    assert by_clade["A|B"].node_age == 0.1
    assert by_clade["A|B"].metadata[0].values_by_taxon == [
        "A=Alpha species",
        "B=Beta species",
    ]
    assert by_clade["A|B"].metadata[1].distinct_values == ["Norway", "Sweden"]

    assert by_clade["A"].node_kind == "tip"
    assert by_clade["A"].support is None
    assert by_clade["A"].node_age == 0.0
    assert by_clade["A"].metadata[0].values_by_taxon == ["A=Alpha species"]


def test_extract_tree_set_clades_reports_one_row_per_tree_clade() -> None:
    report = extract_tree_set_clades(tree_fixture("example_tree_set_left.nwk"))

    assert report.tree_count == 3
    assert len(report.rows) == 21
    assert sorted({row.tree_index for row in report.rows}) == [1, 2, 3]

    ab_rows = [row for row in report.rows if row.clade_id == "A|B"]
    assert len(ab_rows) == 2
    assert [row.tree_index for row in ab_rows] == [1, 2]
    assert all(row.node_kind == "internal" for row in ab_rows)

    ac_rows = [row for row in report.rows if row.clade_id == "A|C"]
    assert len(ac_rows) == 1
    assert ac_rows[0].tree_index == 3


def test_extract_tree_set_clades_requires_identical_taxon_sets() -> None:
    with pytest.raises(InvalidAlignmentError):
        extract_tree_set_clades(tree_fixture("example_tree_set_mismatched.nwk"))
