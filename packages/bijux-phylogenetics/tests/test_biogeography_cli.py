from __future__ import annotations

import json
from pathlib import Path

import pytest

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


def test_biogeography_model_cli_can_export_review(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    nodes_path = tmp_path / "nodes.tsv"
    rates_path = tmp_path / "rates.tsv"
    events_path = tmp_path / "events.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"

    exit_code = main(
        [
            "biogeography",
            "model",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "ard",
            "--summary-out",
            str(summary_path),
            "--nodes-out",
            str(nodes_path),
            "--rates-out",
            str(rates_path),
            "--events-out",
            str(events_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "ard"
    assert payload["metrics"]["transition_rate_row_count"] > 0
    assert payload["metrics"]["changed_branch_count"] >= 0
    assert "root_region" in summary_path.read_text(encoding="utf-8")
    assert "most_likely_region" in nodes_path.read_text(encoding="utf-8")
    assert "source_region" in rates_path.read_text(encoding="utf-8")
    assert "strongly_supported" in events_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_biogeography_model_cli_accepts_region_vocabulary(capsys) -> None:
    exit_code = main(
        [
            "biogeography",
            "model",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "sym",
            "--allowed-regions",
            "north,south,island",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "sym"
    assert payload["metrics"]["observed_region_count"] == 3


@pytest.mark.slow
def test_biogeography_constrained_cli_can_export_review(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    fits_path = tmp_path / "fits.tsv"
    transitions_path = tmp_path / "transitions.tsv"
    unsupported_path = tmp_path / "unsupported.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"

    exit_code = main(
        [
            "biogeography",
            "constrained",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            str(fixture("example_geographic_adjacency.tsv")),
            "--trait",
            "region",
            "--model",
            "ard",
            "--summary-out",
            str(summary_path),
            "--fits-out",
            str(fits_path),
            "--transitions-out",
            str(transitions_path),
            "--unsupported-out",
            str(unsupported_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "ard"
    assert payload["metrics"]["allowed_transition_count"] > 0
    assert payload["metrics"]["forbidden_transition_count"] > 0
    assert "preferred_constraint" in summary_path.read_text(encoding="utf-8")
    assert "constraint_mode" in fits_path.read_text(encoding="utf-8")
    assert "transition_allowed" in transitions_path.read_text(encoding="utf-8")
    assert "claim_resolved" in unsupported_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_biogeography_time_stratified_cli_can_export_review(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    matrix_path = tmp_path / "matrix.tsv"
    branches_path = tmp_path / "branches.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"

    exit_code = main(
        [
            "biogeography",
            "time-stratified",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "ard",
            "--time-bin",
            "early:0.0:0.1",
            "--time-bin",
            "late:0.1:0.3",
            "--summary-out",
            str(summary_path),
            "--matrix-out",
            str(matrix_path),
            "--branches-out",
            str(branches_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "ard"
    assert payload["metrics"]["time_bin_count"] == 2
    assert payload["metrics"]["matrix_row_count"] > 0
    assert "tree_depth" in summary_path.read_text(encoding="utf-8")
    assert "time_stratified_rate" in matrix_path.read_text(encoding="utf-8")
    assert "allocated_transition_weight" in branches_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_biogeography_time_stratified_cli_rejects_invalid_time_bins(capsys) -> None:
    with pytest.raises(SystemExit) as error:
        main(
            [
                "biogeography",
                "time-stratified",
                str(fixture("example_tree.nwk")),
                str(fixture("example_traits_geography.tsv")),
                "--trait",
                "region",
                "--time-bin",
                "broken",
            ]
        )
    error_text = capsys.readouterr().err

    assert error.value.code == 2
    assert "LABEL:START:END" in error_text


def test_biogeography_sampling_bias_cli_can_export_review(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    regions_path = tmp_path / "regions.tsv"
    nodes_path = tmp_path / "nodes.tsv"
    transitions_path = tmp_path / "transitions.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"

    exit_code = main(
        [
            "biogeography",
            "sampling-bias",
            str(fixture("example_tree_six_taxa.nwk")),
            str(fixture("example_traits_geography_biased.tsv")),
            "--trait",
            "region",
            "--model",
            "ard",
            "--weights",
            str(fixture("example_geographic_region_weights.tsv")),
            "--summary-out",
            str(summary_path),
            "--regions-out",
            str(regions_path),
            "--nodes-out",
            str(nodes_path),
            "--transitions-out",
            str(transitions_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "ard"
    assert payload["metrics"]["weighting_mode"] == "explicit"
    assert payload["metrics"]["root_region_changed"] is True
    assert "weighting_mode" in summary_path.read_text(encoding="utf-8")
    assert "weighted_sample_fraction" in regions_path.read_text(encoding="utf-8")
    assert "weighted_region_probabilities" in nodes_path.read_text(encoding="utf-8")
    assert "changed_by_weighting" in transitions_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_biogeography_chronology_cli_can_export_review(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    nodes_path = tmp_path / "nodes.tsv"
    events_path = tmp_path / "events.tsv"
    bins_path = tmp_path / "bins.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"

    exit_code = main(
        [
            "biogeography",
            "chronology",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "ard",
            "--time-bin-count",
            "3",
            "--summary-out",
            str(summary_path),
            "--nodes-out",
            str(nodes_path),
            "--events-out",
            str(events_path),
            "--bins-out",
            str(bins_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["tree_is_time_scaled"] is True
    assert payload["metrics"]["root_age"] == 0.3
    assert payload["metrics"]["time_bin_count"] == 3
    assert "tree_is_time_scaled" in summary_path.read_text(encoding="utf-8")
    assert "age_before_present" in nodes_path.read_text(encoding="utf-8")
    assert "midpoint_age_before_present" in events_path.read_text(encoding="utf-8")
    assert "uncertainty_class" in bins_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_biogeography_events_cli_can_export_single_tree_review(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    events_path = tmp_path / "events.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"

    exit_code = main(
        [
            "biogeography",
            "events",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "ard",
            "--summary-out",
            str(summary_path),
            "--events-out",
            str(events_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["report_mode"] == "single_tree"
    assert payload["metrics"]["event_count"] == 2
    assert "mean_event_support" in summary_path.read_text(encoding="utf-8")
    assert "midpoint_depth" in events_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_biogeography_events_cli_can_export_tree_set_review(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "summary.tsv"
    trees_path = tmp_path / "trees.tsv"
    events_path = tmp_path / "events.tsv"
    event_summaries_path = tmp_path / "event-summaries.tsv"
    exclusions_path = tmp_path / "exclusions.tsv"

    exit_code = main(
        [
            "biogeography",
            "events",
            str(fixture("example_tree_set_left.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            "--model",
            "ard",
            "--tree-set",
            "--summary-out",
            str(summary_path),
            "--trees-out",
            str(trees_path),
            "--events-out",
            str(events_path),
            "--event-summaries-out",
            str(event_summaries_path),
            "--exclusions-out",
            str(exclusions_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["report_mode"] == "tree_set"
    assert payload["metrics"]["kept_tree_count"] == 3
    assert "event_summary_count" in summary_path.read_text(encoding="utf-8")
    assert "rooted_topology_id" in trees_path.read_text(encoding="utf-8")
    assert "midpoint_depth" in events_path.read_text(encoding="utf-8")
    assert "tree_presence_fraction" in event_summaries_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_biogeography_report_cli_can_export_full_review_package(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "biogeography-report"

    exit_code = main(
        [
            "biogeography",
            "report",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits_geography.tsv")),
            "--trait",
            "region",
            str(fixture("example_geographic_region_centroids.tsv")),
            "--model",
            "ard",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["report_kind"] == "biogeography-report-package"
    assert payload["metrics"]["artifact_count"] == 15
    assert payload["metrics"]["event_count"] == 2
    assert payload["metrics"]["publication_ready"] is True
    assert payload["metrics"]["caption_ready"] is True
    assert payload["metrics"]["rendered_internal_pie_count"] > 0
    assert (out_dir / "biogeography-report.html").exists()
    assert (out_dir / "ancestral-region-tree.svg").exists()
    assert (out_dir / "geographic-region-map.html").exists()
    assert (out_dir / "figure-legend.tsv").exists()
    assert (out_dir / "figure-caption.md").exists()
    assert (out_dir / "summary.tsv").exists()
    assert (out_dir / "region-counts.tsv").exists()
    assert (out_dir / "ancestral-regions.tsv").exists()
    assert (out_dir / "transition-matrix.tsv").exists()
    assert (out_dir / "event-table.tsv").exists()
    assert (out_dir / "map-markers.tsv").exists()
    assert (out_dir / "map-lines.tsv").exists()
    assert (out_dir / "exclusions.tsv").exists()
    assert (out_dir / "biogeography-report.manifest.json").exists()
    assert (out_dir / "figure-reproducibility.manifest.json").exists()
