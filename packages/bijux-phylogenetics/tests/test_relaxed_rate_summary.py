from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.dating as dating_api
from bijux_phylogenetics.phylo.dating import (
    summarize_relaxed_rate_branches_from_paths,
)
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "metadata")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_package_dating_gateway_exports_relaxed_rate_summary_surface() -> None:
    assert (
        dating_api.summarize_relaxed_rate_branches_from_paths
        is summarize_relaxed_rate_branches_from_paths
    )


def test_relaxed_rate_summary_reports_branch_rates_and_outlier() -> None:
    report = summarize_relaxed_rate_branches_from_paths(
        fixture("relaxed_rate_summary_substitution_tree_4_taxa.nwk"),
        fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk"),
    )

    assert report.substitution_tree_newick == "(((A:1.2,B:0.2):0.4,C:0.4):0.4,D:0.8);"
    assert report.dated_tree_newick == "(((A:2,B:2):2,C:4):4,D:8);"
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.branch_count == 6
    assert report.substitution_tree_path == str(
        fixture("relaxed_rate_summary_substitution_tree_4_taxa.nwk")
    )
    assert report.dated_tree_path == str(
        fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk")
    )
    assert report.outlier_threshold == pytest.approx(2.0, abs=1e-12)
    assert report.mean_branch_rate == pytest.approx(0.2, abs=1e-12)
    assert report.standard_deviation_branch_rate == pytest.approx(
        0.18257418583505536,
        abs=1e-15,
    )
    assert report.minimum_branch_rate == pytest.approx(0.1, abs=1e-12)
    assert report.maximum_branch_rate == pytest.approx(0.6, abs=1e-12)
    assert report.outlier_count == 1

    rows_by_descendant_taxa = {
        tuple(row.descendant_taxa): row for row in report.branch_rows
    }
    assert rows_by_descendant_taxa[("A",)].branch_rate == pytest.approx(0.6, abs=1e-12)
    assert rows_by_descendant_taxa[("A",)].rate_z_score == pytest.approx(
        2.1908902300206643,
        abs=1e-15,
    )
    assert rows_by_descendant_taxa[("A",)].outlier is True
    assert rows_by_descendant_taxa[("A", "B")].branch_rate == pytest.approx(
        0.2,
        abs=1e-12,
    )
    assert rows_by_descendant_taxa[("A", "B")].rate_z_score == pytest.approx(
        0.0,
        abs=1e-15,
    )
    assert rows_by_descendant_taxa[("D",)].branch_rate == pytest.approx(0.1, abs=1e-12)
    assert rows_by_descendant_taxa[("D",)].outlier is False
    assert len(report.outlier_rows) == 1
    assert report.outlier_rows[0].child_name == "A"
    assert report.outlier_rows[0].descendant_taxa == ["A"]


def test_relaxed_rate_summary_reports_zero_z_scores_on_constant_rate_fixture() -> None:
    report = summarize_relaxed_rate_branches_from_paths(
        fixture("least_squares_dating_substitution_tree_4_taxa.nwk"),
        fixture("least_squares_dating_time_tree_4_taxa.nwk"),
    )

    assert report.mean_branch_rate == pytest.approx(0.25, abs=1e-12)
    assert report.standard_deviation_branch_rate == pytest.approx(0.0, abs=1e-12)
    assert report.outlier_count == 0
    assert all(
        row.branch_rate == pytest.approx(0.25, abs=1e-12) for row in report.branch_rows
    )
    assert all(
        row.rate_z_score == pytest.approx(0.0, abs=1e-12) for row in report.branch_rows
    )
    assert report.outlier_rows == []


def test_relaxed_rate_summary_rejects_mismatched_topology(tmp_path: Path) -> None:
    mismatched_dated_tree_path = tmp_path / "mismatched-dated-tree.nwk"
    mismatched_dated_tree_path.write_text(
        "(((A:2,C:2):2,B:4):4,D:8);",
        encoding="utf-8",
    )

    with pytest.raises(
        PhylogeneticsError,
        match="requires identical rooted topology",
    ):
        summarize_relaxed_rate_branches_from_paths(
            fixture("relaxed_rate_summary_substitution_tree_4_taxa.nwk"),
            mismatched_dated_tree_path,
        )


def test_relaxed_rate_summary_rejects_nonpositive_dated_durations(
    tmp_path: Path,
) -> None:
    zero_duration_tree_path = tmp_path / "zero-duration-tree.nwk"
    zero_duration_tree_path.write_text(
        "(((A:0,B:2):2,C:4):4,D:8);",
        encoding="utf-8",
    )

    with pytest.raises(
        InvalidBranchLengthError,
        match="requires strictly positive dated branch durations",
    ):
        summarize_relaxed_rate_branches_from_paths(
            fixture("relaxed_rate_summary_substitution_tree_4_taxa.nwk"),
            zero_duration_tree_path,
        )
