from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.trees import (
    TREE_SET_SPLIT_FREQUENCY_POLICIES,
    compute_tree_set_split_frequency_table,
    write_tree_set_split_frequency_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_compute_tree_set_split_frequency_table_counts_rooted_signatures() -> None:
    report = compute_tree_set_split_frequency_table(
        fixture("example_tree_set_rooting_only_difference.nwk"),
        split_policy="rooted",
    )

    assert report.split_policy == "rooted"
    assert report.tree_count == 2
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert [
        (row.split, row.tree_count, row.frequency) for row in report.split_frequencies
    ] == [
        ("A|B", 2, 1.0),
        ("A|B|C", 1, 0.5),
        ("C|D", 1, 0.5),
    ]


def test_compute_tree_set_split_frequency_table_counts_unrooted_signatures() -> None:
    report = compute_tree_set_split_frequency_table(
        fixture("example_tree_set_rooting_only_difference.nwk"),
        split_policy="unrooted",
    )

    assert report.split_policy == "unrooted"
    assert report.tree_count == 2
    assert [
        (row.split, row.tree_count, row.frequency) for row in report.split_frequencies
    ] == [("A|B", 2, 1.0)]


def test_compute_tree_set_split_frequency_table_distinguishes_rooting_policies() -> (
    None
):
    rooted = compute_tree_set_split_frequency_table(
        fixture("example_tree_set_rooting_only_difference.nwk"),
        split_policy="rooted",
    )
    unrooted = compute_tree_set_split_frequency_table(
        fixture("example_tree_set_rooting_only_difference.nwk"),
        split_policy="unrooted",
    )

    assert TREE_SET_SPLIT_FREQUENCY_POLICIES == ("rooted", "unrooted")
    assert [row.split for row in rooted.split_frequencies] == ["A|B", "A|B|C", "C|D"]
    assert [row.split for row in unrooted.split_frequencies] == ["A|B"]


def test_write_tree_set_split_frequency_table_writes_expected_columns(
    tmp_path: Path,
) -> None:
    report = compute_tree_set_split_frequency_table(
        fixture("example_tree_set_rooting_only_difference.nwk"),
        split_policy="rooted",
    )

    output_path = tmp_path / "tree-set-split-frequencies.tsv"
    write_tree_set_split_frequency_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()

    assert lines[0] == "split\ttree_count\tfrequency\tsplit_policy"
    assert "A|B\t2\t1\trooted" in lines[1:]


def test_compute_tree_set_split_frequency_table_rejects_unknown_policy() -> None:
    with pytest.raises(ValueError, match="split_policy must be one of"):
        compute_tree_set_split_frequency_table(
            fixture("example_tree_set_rooting_only_difference.nwk"),
            split_policy="tips-only",
        )
