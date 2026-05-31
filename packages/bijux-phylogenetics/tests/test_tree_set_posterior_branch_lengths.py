from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.trees import (
    summarize_posterior_branch_lengths,
    write_posterior_branch_length_summary_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_summarize_posterior_branch_lengths_matches_clades_across_topologies() -> None:
    report = summarize_posterior_branch_lengths(
        fixture("posterior_branch_length_summary_tree_set.nwk")
    )
    rows = {row.clade: row for row in report.rows}

    assert report.tree_count == 4
    assert report.shared_taxa == ["A", "B", "C", "D", "E"]
    assert report.hpd_mass == 0.95
    assert [row.clade for row in report.rows[:3]] == ["A|B", "D|E", "C|D|E"]

    a_b = rows["A|B"]
    assert a_b.matched_tree_count == 4
    assert a_b.posterior_tree_count == 4
    assert a_b.clade_frequency == 1.0
    assert a_b.mean_branch_length == pytest.approx(0.25)
    assert a_b.median_branch_length == pytest.approx(0.25)
    assert a_b.hpd_95_lower == pytest.approx(0.1)
    assert a_b.hpd_95_upper == pytest.approx(0.4)
    assert a_b.effective_sample_size is not None
    assert 1.0 <= a_b.effective_sample_size <= 4.0

    d_e = rows["D|E"]
    assert d_e.matched_tree_count == 3
    assert d_e.clade_frequency == pytest.approx(0.75)
    assert d_e.mean_branch_length == pytest.approx(0.5)
    assert d_e.median_branch_length == pytest.approx(0.5)
    assert d_e.hpd_95_lower == pytest.approx(0.4)
    assert d_e.hpd_95_upper == pytest.approx(0.6)

    c_d_e = rows["C|D|E"]
    assert c_d_e.matched_tree_count == 2
    assert c_d_e.clade_frequency == pytest.approx(0.5)
    assert c_d_e.mean_branch_length == pytest.approx(0.75)
    assert c_d_e.median_branch_length == pytest.approx(0.75)


def test_write_posterior_branch_length_summary_table_writes_expected_columns(
    tmp_path: Path,
) -> None:
    report = summarize_posterior_branch_lengths(
        fixture("posterior_branch_length_summary_tree_set.nwk")
    )

    output_path = tmp_path / "posterior-branch-lengths.tsv"
    write_posterior_branch_length_summary_table(output_path, report)

    assert output_path.read_text(encoding="utf-8").splitlines()[:3] == [
        "clade\tmatched_tree_count\tposterior_tree_count\tclade_frequency\tmean_branch_length\tmedian_branch_length\thpd_95_lower\thpd_95_upper\teffective_sample_size",
        "A|B\t4\t4\t1\t0.25\t0.25\t0.1\t0.4\t2.4",
        "D|E\t3\t4\t0.75\t0.5\t0.5\t0.4\t0.6\t3",
    ]


def test_summarize_posterior_branch_lengths_requires_exact_taxa() -> None:
    with pytest.raises(InvalidAlignmentError, match="exact same taxon set"):
        summarize_posterior_branch_lengths(fixture("example_tree_set_mismatched.nwk"))
