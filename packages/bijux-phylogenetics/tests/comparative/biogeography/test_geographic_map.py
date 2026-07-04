from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylogeography import (
    render_geographic_map_html,
    summarize_continuous_phylogeography_map,
    summarize_discrete_region_map,
    write_geographic_map_exclusion_table,
    write_geographic_map_line_table,
    write_geographic_map_marker_table,
    write_geographic_map_summary_table,
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


def test_summarize_continuous_phylogeography_map_supports_depth_filter_and_html(
    tmp_path: Path,
) -> None:
    report = summarize_continuous_phylogeography_map(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogeography.tsv"),
        latitude_column="latitude",
        longitude_column="longitude",
        model="brownian",
        minimum_midpoint_depth=2.0,
    )

    assert report.mode == "continuous"
    assert report.summary.line_count == len(report.line_rows)
    assert 0 < report.summary.visible_line_count < report.summary.line_count
    assert report.summary.time_filter_applied is True
    assert any(row.is_tip for row in report.marker_rows)
    assert any(not row.is_tip for row in report.marker_rows)

    html_path = render_geographic_map_html(
        report,
        out_path=tmp_path / "continuous-map.html",
    ).output_path
    summary_path = write_geographic_map_summary_table(tmp_path / "summary.tsv", report)
    markers_path = write_geographic_map_marker_table(tmp_path / "markers.tsv", report)
    lines_path = write_geographic_map_line_table(tmp_path / "lines.tsv", report)
    exclusions_path = write_geographic_map_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )

    assert "Continuous Geographic Map Review" in html_path.read_text(encoding="utf-8")
    assert "visible_line_count" in summary_path.read_text(encoding="utf-8")
    assert "marker_kind" in markers_path.read_text(encoding="utf-8")
    assert "state_transition" in lines_path.read_text(encoding="utf-8")
    assert "subject_kind" in exclusions_path.read_text(encoding="utf-8")


def test_summarize_discrete_region_map_supports_event_filter_and_html(
    tmp_path: Path,
) -> None:
    report = summarize_discrete_region_map(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        centroids_path=fixture("example_geographic_region_centroids.tsv"),
        model="ard",
        maximum_midpoint_depth=0.15,
    )

    assert report.mode == "regions"
    assert report.summary.line_count == 2
    assert report.summary.visible_line_count == 1
    assert any(row.state_label for row in report.marker_rows)
    assert any(row.state_transition == "north->south" for row in report.line_rows)

    html_path = render_geographic_map_html(
        report,
        out_path=tmp_path / "regions-map.html",
    ).output_path

    html_text = html_path.read_text(encoding="utf-8")
    assert "Regional Transition Map Review" in html_text
    assert "north-&gt;south" in html_text or "north->south" in html_text
    assert "midpoint-depth filter" in html_text


def test_summarize_discrete_region_map_tracks_missing_centroids(tmp_path: Path) -> None:
    centroid_path = tmp_path / "centroids.tsv"
    centroid_path.write_text(
        "region\tlatitude\tlongitude\nnorth\t59.0\t18.0\nsouth\t51.0\t10.0\n",
        encoding="utf-8",
    )

    report = summarize_discrete_region_map(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        centroids_path=centroid_path,
        model="er",
    )

    assert report.summary.excluded_record_count >= 1
    assert any(row.reason == "missing-region-centroid" for row in report.exclusion_rows)
