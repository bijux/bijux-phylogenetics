from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.cli import main

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
    assert "tree_presence_fraction" in event_summaries_path.read_text(
        encoding="utf-8"
    )
    assert "reason" in exclusions_path.read_text(encoding="utf-8")
