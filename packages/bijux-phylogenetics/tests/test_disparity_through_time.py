from __future__ import annotations

import csv
import math
from pathlib import Path

from bijux_phylogenetics.comparative.disparity import (
    render_disparity_through_time_svg,
    summarize_continuous_clade_disparity,
    summarize_disparity_through_time,
    write_continuous_clade_disparity_table,
    write_disparity_through_time_bin_table,
    write_disparity_through_time_curve_table,
    write_disparity_through_time_exclusion_table,
    write_disparity_through_time_summary_table,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")

_GEIGER_DTT_TIMES = [
    0.0,
    0.0,
    0.2031936,
    0.2085324,
    0.3994185,
    0.7157151,
    0.7345286,
    0.8212913,
    0.8285889,
    0.8542111,
    0.8616957,
    0.8650965,
    0.9015125,
    0.9060953,
    0.9220957,
    0.9382944,
    0.9687149,
    0.9691681,
    0.9765662,
    0.9838479,
    0.9947175,
    0.9947625,
    0.9973615,
    0.9991831,
]

_GEIGER_DTT_UNIVARIATE = [
    1.0,
    0.512322522,
    0.663785596,
    0.548801416,
    0.553639756,
    0.390989837,
    0.374633159,
    0.506858092,
    0.439899252,
    0.366992203,
    0.251160342,
    0.211947133,
    0.124119626,
    0.123827183,
    0.090213561,
    0.084967964,
    0.084252832,
    0.101002688,
    0.077085679,
    0.046832117,
    0.036559716,
    0.007257707,
    0.013896031,
    0.0,
]

_GEIGER_DTT_MULTIVARIATE = [
    1.0,
    0.507465918,
    0.777560903,
    0.440187451,
    0.442661666,
    0.312776058,
    0.299716585,
    0.40499346,
    0.351331319,
    0.29312038,
    0.200713999,
    0.169413498,
    0.099178564,
    0.09890281,
    0.072088886,
    0.067893208,
    0.067299035,
    0.080644235,
    0.061560397,
    0.037397332,
    0.02918871,
    0.005808167,
    0.011089956,
    0.0,
]


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_summarize_continuous_clade_disparity_reports_geiger_style_root_row() -> None:
    report = summarize_continuous_clade_disparity(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_geiger_continuous_model_panel_twenty_four_taxa.tsv"),
        trait_columns=["ou_truth"],
    )

    assert report.trait_columns == ["ou_truth"]
    assert report.analyzed_taxon_count == 24
    assert len(report.clade_rows) == 23
    assert report.clade_rows[0].ape_node_id == 25
    assert report.clade_rows[0].node_kind == "root"
    assert report.clade_rows[0].descendant_taxon_count == 24
    assert math.isclose(report.root_age, 1.6918111678463008)
    assert math.isclose(
        report.clade_rows[0].disparity, 0.015837416289029, rel_tol=1e-12
    )


def test_summarize_disparity_through_time_matches_geiger_univariate_curve() -> None:
    report = summarize_disparity_through_time(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_geiger_continuous_model_panel_twenty_four_taxa.tsv"),
        trait_columns=["ou_truth"],
    )

    assert report.relative_scaling_applied is True
    assert len(report.curve_rows) == len(_GEIGER_DTT_TIMES)
    for row, expected_time, expected_disparity in zip(
        report.curve_rows,
        _GEIGER_DTT_TIMES,
        _GEIGER_DTT_UNIVARIATE,
        strict=True,
    ):
        assert math.isclose(
            row.relative_time, expected_time, rel_tol=1e-6, abs_tol=1e-6
        )
        assert math.isclose(
            row.relative_disparity,
            expected_disparity,
            rel_tol=1e-9,
            abs_tol=1e-9,
        )


def test_summarize_disparity_through_time_matches_geiger_multivariate_curve() -> None:
    report = summarize_disparity_through_time(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_geiger_continuous_model_panel_twenty_four_taxa.tsv"),
        trait_columns=["ou_truth", "early_burst_truth"],
    )

    assert report.trait_columns == ["ou_truth", "early_burst_truth"]
    assert len(report.curve_rows) == len(_GEIGER_DTT_TIMES)
    for row, expected_time, expected_disparity in zip(
        report.curve_rows,
        _GEIGER_DTT_TIMES,
        _GEIGER_DTT_MULTIVARIATE,
        strict=True,
    ):
        assert math.isclose(
            row.relative_time, expected_time, rel_tol=1e-6, abs_tol=1e-6
        )
        assert math.isclose(
            row.relative_disparity,
            expected_disparity,
            rel_tol=1e-9,
            abs_tol=1e-9,
        )


def test_summarize_disparity_through_time_prunes_missing_traits_and_supports_bins(
    tmp_path: Path,
) -> None:
    table_path = tmp_path / "dtt-missing.tsv"
    table_path.write_text(
        fixture("example_traits_geiger_continuous_model_panel_twenty_four_taxa.tsv")
        .read_text(encoding="utf-8")
        .replace("Phy10\t0.631983", "Phy10\t"),
        encoding="utf-8",
    )

    report = summarize_disparity_through_time(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        table_path,
        trait_columns=["ou_truth", "early_burst_truth"],
        time_bin_count=4,
    )

    assert report.analyzed_taxon_count == 23
    assert any(
        row.taxon == "Phy10" and row.reason == "missing_trait_value"
        for row in report.excluded_taxa
    )
    assert len(report.time_bin_rows) == 4
    assert sum(row.point_count for row in report.time_bin_rows) == len(
        report.curve_rows
    )


def test_disparity_through_time_writers_emit_tables_and_svg(tmp_path: Path) -> None:
    report = summarize_disparity_through_time(
        fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("example_traits_geiger_continuous_model_panel_twenty_four_taxa.tsv"),
        trait_columns=["ou_truth", "early_burst_truth"],
        time_bin_count=3,
    )
    summary_path = tmp_path / "dtt-summary.tsv"
    curve_path = tmp_path / "dtt-curve.tsv"
    clade_path = tmp_path / "clade-disparity.tsv"
    exclusion_path = tmp_path / "dtt-excluded.tsv"
    bin_path = tmp_path / "dtt-bins.tsv"
    svg_path = tmp_path / "dtt.svg"

    write_disparity_through_time_summary_table(summary_path, report)
    write_disparity_through_time_curve_table(curve_path, report)
    write_continuous_clade_disparity_table(
        clade_path,
        summarize_continuous_clade_disparity(
            fixture("example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
            fixture(
                "example_traits_geiger_continuous_model_panel_twenty_four_taxa.tsv"
            ),
            trait_columns=["ou_truth", "early_burst_truth"],
        ),
    )
    write_disparity_through_time_exclusion_table(exclusion_path, report)
    write_disparity_through_time_bin_table(bin_path, report)
    plotted = render_disparity_through_time_svg(svg_path, report)

    assert plotted == len(report.curve_rows)
    assert "<svg" in svg_path.read_text(encoding="utf-8")
    with curve_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == len(report.curve_rows)
    assert summary_path.exists()
    assert clade_path.exists()
    assert exclusion_path.exists()
    assert bin_path.exists()
