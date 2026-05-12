from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylogeography import (
    render_coordinate_movement_visualization,
    summarize_continuous_phylogeography,
    write_coordinate_estimate_table,
    write_coordinate_movement_branch_table,
    write_coordinate_movement_exclusion_table,
    write_coordinate_movement_outlier_table,
    write_coordinate_movement_summary_table,
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


def test_summarize_continuous_phylogeography_supports_brownian_and_ou() -> None:
    for model in ("brownian", "ou"):
        report = summarize_continuous_phylogeography(
            fixture("example_tree_six_taxa.nwk"),
            fixture("example_traits_phylogeography.tsv"),
            latitude_column="latitude",
            longitude_column="longitude",
            model=model,
            alpha=0.5,
        )

        assert report.model == model
        assert report.summary.model == model
        assert report.estimate_rows
        assert report.branch_rows
        assert report.outlier_rows


def test_summarize_continuous_phylogeography_reports_uncertainty_and_outliers() -> None:
    report = summarize_continuous_phylogeography(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogeography.tsv"),
        latitude_column="latitude",
        longitude_column="longitude",
        model="brownian",
    )

    root_row = next(row for row in report.estimate_rows if row.is_root)

    assert report.summary.internal_node_count == 5
    assert report.summary.flagged_branch_count >= 1
    assert root_row.radial_standard_error_km >= 0.0
    assert any(row.outlier_jump for row in report.branch_rows)


def test_summarize_continuous_phylogeography_tracks_coordinate_exclusions(
    tmp_path: Path,
) -> None:
    table_path = tmp_path / "coordinates.tsv"
    table_path.write_text(
        (
            "taxon\tlatitude\tlongitude\n"
            "A\t59.0\t18.0\n"
            "B\t\t19.0\n"
            "C\t95.0\t13.0\n"
            "D\t48.0\tnot-a-number\n"
            "E\t-33.0\t151.0\n"
            "Ghost\t10.0\t10.0\n"
        ),
        encoding="utf-8",
    )

    report = summarize_continuous_phylogeography(
        fixture("example_tree_six_taxa.nwk"),
        table_path,
        latitude_column="latitude",
        longitude_column="longitude",
        model="brownian",
    )

    assert {(row.taxon, row.reason) for row in report.exclusion_rows} == {
        ("B", "missing-coordinate"),
        ("C", "latitude-out-of-range"),
        ("D", "non-numeric-coordinate"),
        ("Ghost", "taxon-not-in-tree"),
    }


def test_render_coordinate_movement_visualization_and_ledgers(
    tmp_path: Path,
) -> None:
    report = summarize_continuous_phylogeography(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogeography.tsv"),
        latitude_column="latitude",
        longitude_column="longitude",
        model="brownian",
    )

    summary_path = write_coordinate_movement_summary_table(tmp_path / "summary.tsv", report)
    estimates_path = write_coordinate_estimate_table(tmp_path / "estimates.tsv", report)
    branches_path = write_coordinate_movement_branch_table(tmp_path / "branches.tsv", report)
    outliers_path = write_coordinate_movement_outlier_table(tmp_path / "outliers.tsv", report)
    exclusions_path = write_coordinate_movement_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )
    svg_path = tmp_path / "movement.svg"
    html_path = tmp_path / "movement.html"
    svg_result = render_coordinate_movement_visualization(report, out_path=svg_path)
    html_result = render_coordinate_movement_visualization(report, out_path=html_path)

    assert "root_latitude" in summary_path.read_text(encoding="utf-8")
    assert "radial_standard_error_km" in estimates_path.read_text(encoding="utf-8")
    assert "great_circle_km" in branches_path.read_text(encoding="utf-8")
    assert "flag_codes" in outliers_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")
    assert svg_result.output_path == svg_path
    assert html_result.output_path == html_path
    assert "<svg" in svg_path.read_text(encoding="utf-8")
    assert "Coordinate-Space Movement Review" in html_path.read_text(encoding="utf-8")
