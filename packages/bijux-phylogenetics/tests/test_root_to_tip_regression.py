from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.diagnostics.root_to_tip import (
    diagnose_root_to_tip_regression,
)
from bijux_phylogenetics.runtime.errors import (
    MetadataJoinError,
    PhylogeneticsError,
    UnrootedTreeError,
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


def test_diagnose_root_to_tip_regression_reports_fit_and_outlier() -> None:
    report = diagnose_root_to_tip_regression(
        fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk"),
        fixture("root_to_tip_regression_dates_7_taxa.tsv"),
    )

    assert report.source_format == "newick"
    assert report.taxon_column == "taxon"
    assert report.date_column == "date"
    assert report.tip_count == 7
    assert report.slope == pytest.approx(2.392857142857143, abs=1e-12)
    assert report.intercept == pytest.approx(-1.821428571428572, abs=1e-12)
    assert report.r_squared == pytest.approx(0.6390945330296127, abs=1e-12)
    assert report.residual_mean_square == pytest.approx(18.107142857142858, abs=1e-12)
    assert [row.tip for row in report.rows] == ["A", "B", "C", "D", "E", "F", "G"]
    assert report.rows[-1].outlier is True
    assert report.rows[-1].studentized_residual == pytest.approx(
        2.23606797749979,
        abs=1e-12,
    )
    assert [row.tip for row in report.outliers] == ["G"]
    assert report.outliers[0].residual == pytest.approx(6.964285714285715, abs=1e-12)


def test_diagnose_root_to_tip_regression_rejects_unrooted_trees(tmp_path: Path) -> None:
    metadata_path = tmp_path / "dates.tsv"
    metadata_path.write_text(
        "taxon\tdate\nA\t0\nB\t1\nC\t2\nD\t3\n",
        encoding="utf-8",
    )

    with pytest.raises(UnrootedTreeError, match="tree is not rooted"):
        diagnose_root_to_tip_regression(
            fixture("example_tree_unrooted.nwk"),
            metadata_path,
        )


def test_diagnose_root_to_tip_regression_rejects_missing_tree_taxa(
    tmp_path: Path,
) -> None:
    metadata_path = tmp_path / "incomplete_dates.tsv"
    metadata_path.write_text(
        "taxon\tdate\nA\t0\nB\t1\nC\t2\nD\t3\nE\t4\nF\t5\n",
        encoding="utf-8",
    )

    with pytest.raises(
        MetadataJoinError,
        match="sampling-time table is missing tree taxa: G",
    ):
        diagnose_root_to_tip_regression(
            fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk"),
            metadata_path,
        )


def test_diagnose_root_to_tip_regression_requires_variation_in_sampling_time(
    tmp_path: Path,
) -> None:
    metadata_path = tmp_path / "constant_dates.tsv"
    metadata_path.write_text(
        ("taxon\tdate\nA\t4\nB\t4\nC\t4\nD\t4\nE\t4\nF\t4\nG\t4\n"),
        encoding="utf-8",
    )

    with pytest.raises(
        PhylogeneticsError,
        match="root-to-tip regression requires variation in sampling time",
    ):
        diagnose_root_to_tip_regression(
            fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk"),
            metadata_path,
        )
