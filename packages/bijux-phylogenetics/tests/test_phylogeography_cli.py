from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

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


def test_phylogeography_coordinates_cli_can_export_review(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    estimates_path = tmp_path / "estimates.tsv"
    branches_path = tmp_path / "branches.tsv"
    outliers_path = tmp_path / "outliers.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"
    visualization_path = tmp_path / "movement.svg"

    exit_code = main(
        [
            "phylogeography",
            "coordinates",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_phylogeography.tsv")),
            "--latitude-column",
            "latitude",
            "--longitude-column",
            "longitude",
            "--model",
            "brownian",
            "--summary-out",
            str(summary_path),
            "--estimates-out",
            str(estimates_path),
            "--branches-out",
            str(branches_path),
            "--outliers-out",
            str(outliers_path),
            "--exclusions-out",
            str(exclusions_path),
            "--visualization-out",
            str(visualization_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "brownian"
    assert payload["metrics"]["flagged_branch_count"] >= 1
    assert "root_latitude" in summary_path.read_text(encoding="utf-8")
    assert "radial_standard_error_km" in estimates_path.read_text(encoding="utf-8")
    assert "great_circle_km" in branches_path.read_text(encoding="utf-8")
    assert "flag_codes" in outliers_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")
    assert "<svg" in visualization_path.read_text(encoding="utf-8")


def test_phylogeography_coordinates_map_cli_can_export_html_map(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    markers_path = tmp_path / "markers.tsv"
    lines_path = tmp_path / "lines.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"
    html_path = tmp_path / "map.html"

    exit_code = main(
        [
            "phylogeography",
            "coordinates-map",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_phylogeography.tsv")),
            "--latitude-column",
            "latitude",
            "--longitude-column",
            "longitude",
            "--minimum-midpoint-depth",
            "2.0",
            "--summary-out",
            str(summary_path),
            "--markers-out",
            str(markers_path),
            "--lines-out",
            str(lines_path),
            "--exclusions-out",
            str(exclusions_path),
            "--html-out",
            str(html_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["map_mode"] == "continuous"
    assert payload["metrics"]["time_filter_applied"] is True
    assert payload["metrics"]["visible_line_count"] < payload["metrics"]["line_count"]
    assert "visible_line_count" in summary_path.read_text(encoding="utf-8")
    assert "marker_kind" in markers_path.read_text(encoding="utf-8")
    assert "state_transition" in lines_path.read_text(encoding="utf-8")
    assert "subject_kind" in exclusions_path.read_text(encoding="utf-8")
    assert "Continuous Geographic Map Review" in html_path.read_text(encoding="utf-8")


def test_phylogeography_regions_map_cli_can_export_html_map(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    markers_path = tmp_path / "markers.tsv"
    lines_path = tmp_path / "lines.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"
    html_path = tmp_path / "map.html"

    exit_code = main(
        [
            "phylogeography",
            "regions-map",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--centroids",
            str(fixture("example_geographic_region_centroids.tsv")),
            "--model",
            "ard",
            "--maximum-midpoint-depth",
            "0.15",
            "--summary-out",
            str(summary_path),
            "--markers-out",
            str(markers_path),
            "--lines-out",
            str(lines_path),
            "--exclusions-out",
            str(exclusions_path),
            "--html-out",
            str(html_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["map_mode"] == "regions"
    assert payload["metrics"]["model"] == "ard"
    assert payload["metrics"]["visible_line_count"] == 1
    assert "time_filter_applied" in summary_path.read_text(encoding="utf-8")
    assert "state_label" in markers_path.read_text(encoding="utf-8")
    assert "north->south" in lines_path.read_text(encoding="utf-8")
    assert "subject_kind" in exclusions_path.read_text(encoding="utf-8")
    assert "Regional Transition Map Review" in html_path.read_text(encoding="utf-8")
