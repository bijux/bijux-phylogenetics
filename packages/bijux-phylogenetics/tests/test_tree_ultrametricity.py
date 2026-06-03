from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
    assess_tree_ultrametricity,
    write_tree_ultrametric_table,
)

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


def test_assess_tree_ultrametricity_matches_exact_ultrametric_tree() -> None:
    report = assess_tree_ultrametricity(fixture("example_tree.nwk"))

    assert report.tip_labels == ["A", "B", "C", "D"]
    assert report.rooted is True
    assert report.option == 1
    assert report.criterion_name == "scaled-range"
    assert report.ultrametric is True
    assert report.criterion_value == 0.0
    assert report.max_tip_depth_deviation == 0.0
    assert report.offending_taxa == []
    assert report.rows[0].node_id == 1
    assert math.isclose(report.rows[0].root_to_tip_depth, 0.3, abs_tol=1e-12)


def test_assess_tree_ultrametricity_accepts_near_ultrametric_tree_with_ape_tolerance() -> (
    None
):
    report = assess_tree_ultrametricity(fixture("example_tree_near_ultrametric.nwk"))

    assert report.ultrametric is True
    assert report.tolerance == APE_ULTRAMETRIC_TOLERANCE
    assert report.max_tip_depth_deviation == pytest.approx(1e-9, abs=1e-15)
    assert report.criterion_value < report.tolerance
    assert report.offending_taxa == ["A", "B", "C", "D"]


def test_assess_tree_ultrametricity_rejects_near_ultrametric_tree_under_tight_tolerance() -> (
    None
):
    report = assess_tree_ultrametricity(
        fixture("example_tree_near_ultrametric.nwk"),
        tolerance=1e-12,
    )

    assert report.ultrametric is False
    assert report.criterion_value > report.tolerance
    assert report.offending_taxa == ["A", "B", "C", "D"]


def test_assess_tree_ultrametricity_reports_non_ultrametric_extremes() -> None:
    report = assess_tree_ultrametricity(fixture("example_tree_ladderized.nwk"))

    assert report.ultrametric is False
    assert report.max_tip_depth_deviation == pytest.approx(0.2, abs=1e-12)
    assert report.offending_taxa == ["A", "B", "D"]
    assert report.root_age == pytest.approx(0.3, abs=1e-12)


def test_assess_tree_ultrametricity_supports_variance_option() -> None:
    report = assess_tree_ultrametricity(
        fixture("example_tree.nwk"),
        option=2,
    )

    assert report.option == 2
    assert report.criterion_name == "variance"
    assert report.criterion_value == 0.0
    assert report.ultrametric is True


def test_assess_tree_ultrametricity_rejects_unknown_option() -> None:
    with pytest.raises(ValueError, match="option must be 1 or 2"):
        assess_tree_ultrametricity(fixture("example_tree.nwk"), option=3)


def test_write_tree_ultrametric_table_preserves_tip_order(tmp_path: Path) -> None:
    report = assess_tree_ultrametricity(fixture("example_tree_near_ultrametric.nwk"))
    output_path = tmp_path / "ultrametric-diagnostics.tsv"

    write_tree_ultrametric_table(output_path, report)

    rows = output_path.read_text(encoding="utf-8").splitlines()
    assert rows[0].startswith("node_id\ttip_label\troot_to_tip_depth")
    assert rows[1].startswith("1\tA\t0.300000001")
    assert rows[2].startswith("2\tB\t0.3")
