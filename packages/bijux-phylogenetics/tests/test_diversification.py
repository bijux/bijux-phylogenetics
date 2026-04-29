from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.diversification import (
    compare_diversification_models,
    compute_lineage_through_time_curve,
    detect_incomplete_taxon_sampling_metadata,
    detect_diversification_outlier_clades,
    estimate_diversification_rate,
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


def test_detect_incomplete_taxon_sampling_metadata_reports_missing_and_invalid_rows() -> None:
    report = detect_incomplete_taxon_sampling_metadata(
        fixture("example_tree.nwk"),
        fixture("example_sampling_fractions_incomplete.tsv"),
    )

    assert report.complete is False
    assert report.sampling_column == "sampling_fraction"
    assert report.missing_taxa == ["D"]
    assert [issue.code for issue in report.invalid_rows] == [
        "missing-sampling-fraction",
        "out-of-range-sampling-fraction",
    ]


def test_estimate_diversification_rate_applies_sampling_correction() -> None:
    report = estimate_diversification_rate(
        fixture("example_tree.nwk"),
        metadata_path=fixture("example_sampling_fractions.tsv"),
        model="birth-death",
    )

    assert report.model == "birth-death"
    assert report.crown_age == 0.3
    assert report.observed_tip_count == 4
    assert report.sampling_fraction == 0.75
    assert report.corrected_tip_count == 5.33333333333333
    assert report.birth_rate >= report.net_diversification_rate
    assert report.aic > 0.0


def test_compare_diversification_models_returns_aic_ranked_rows() -> None:
    report = compare_diversification_models(
        fixture("example_tree.nwk"),
        metadata_path=fixture("example_sampling_fractions.tsv"),
    )

    assert report.better_model in {"yule", "birth-death"}
    assert [row.model for row in report.rows] == ["yule", "birth-death"]
    assert all(row.aic > 0.0 for row in report.rows)


def test_detect_diversification_outlier_clades_flags_high_and_low_clades() -> None:
    report = detect_diversification_outlier_clades(fixture("example_tree.nwk"))

    assert report.global_rate > 0.0
    assert {row.classification for row in report.observations} == {"baseline", "high", "low"}
    assert [row.node for row in report.high_diversification_clades] == ["A|B"]
    assert [row.node for row in report.low_diversification_clades] == ["C|D"]
