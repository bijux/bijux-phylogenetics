from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.diversification import (
    compute_lineage_through_time_curve,
    inspect_diversification_time_tree,
    validate_time_tree_for_diversification,
    write_lineage_through_time_table,
)
from bijux_phylogenetics.errors import DiversificationAnalysisError, UnrootedTreeError


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


def test_validate_time_tree_for_diversification_reports_root_age() -> None:
    report = validate_time_tree_for_diversification(fixture("example_tree.nwk"))

    assert report.rooted is True
    assert report.ultrametric is True
    assert report.branch_length_status == "complete"
    assert report.tip_count == 4
    assert report.root_age == 0.3


def test_compute_lineage_through_time_curve_tracks_lineage_increases() -> None:
    report = compute_lineage_through_time_curve(fixture("example_tree.nwk"))

    assert [(point.time_before_present, point.lineage_count) for point in report.points] == [
        (0.3, 2),
        (0.2, 3),
        (0.1, 4),
        (0.0, 4),
    ]


def test_write_lineage_through_time_table_exports_curve(tmp_path: Path) -> None:
    output_path = tmp_path / "ltt.tsv"
    report = compute_lineage_through_time_curve(fixture("example_tree.nwk"))

    write_lineage_through_time_table(output_path, report)

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "node\ttime_before_present\tlineage_count\tevent",
        "A|B|C|D\t0.3\t2\troot",
        "C|D\t0.2\t3\tspeciation",
        "A|B\t0.1\t4\tspeciation",
        "present\t0\t4\tpresent",
    ]


def test_inspect_diversification_time_tree_rejects_invalid_time_tree() -> None:
    with pytest.raises(DiversificationAnalysisError):
        inspect_diversification_time_tree(fixture("example_tree_no_lengths.nwk"))

    with pytest.raises(UnrootedTreeError):
        validate_time_tree_for_diversification(fixture("example_tree_unrooted.nwk"))
