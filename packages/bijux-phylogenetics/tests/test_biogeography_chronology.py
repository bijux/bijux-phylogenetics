from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.biogeography.migration import (
    summarize_biogeographic_transition_chronology,
    write_dated_biogeography_event_table,
    write_dated_biogeography_exclusion_table,
    write_dated_biogeography_node_table,
    write_dated_biogeography_summary_table,
    write_dated_biogeography_time_bin_table,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

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


def test_summarize_biogeographic_transition_chronology_reports_node_ages_and_bins() -> (
    None
):
    report = summarize_biogeographic_transition_chronology(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="ard",
        time_bin_count=3,
    )

    assert report.summary.tree_is_time_scaled is True
    assert report.summary.root_age == 0.3
    assert report.summary.event_count == 2
    assert report.summary.time_bin_count == 3
    assert report.summary.empty_time_bin_count == 1
    assert len(report.node_rows) == 7
    assert any(
        row.is_root and row.age_before_present == 0.3 for row in report.node_rows
    )
    assert any(row.is_tip and row.age_before_present == 0.0 for row in report.node_rows)
    assert {row.time_bin_label for row in report.event_rows} == {"0.1-0.2", "0.2-0.3"}
    assert any(row.uncertainty_class == "no_events" for row in report.time_bin_rows)


def test_summarize_biogeographic_transition_chronology_rejects_non_ultrametric_tree() -> (
    None
):
    with pytest.raises(AncestralReconstructionError) as error:
        summarize_biogeographic_transition_chronology(
            fixture("example_tree_long_branch.nwk"),
            fixture("example_traits_geography.tsv"),
            trait="region",
            model="er",
        )

    assert "ultrametric time tree" in str(error.value)


def test_write_dated_biogeography_tables_emit_expected_ledgers(tmp_path: Path) -> None:
    report = summarize_biogeographic_transition_chronology(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="sym",
        time_bin_count=3,
    )

    summary_path = write_dated_biogeography_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    nodes_path = write_dated_biogeography_node_table(
        tmp_path / "nodes.tsv",
        report,
    )
    events_path = write_dated_biogeography_event_table(
        tmp_path / "events.tsv",
        report,
    )
    bins_path = write_dated_biogeography_time_bin_table(
        tmp_path / "bins.tsv",
        report,
    )
    exclusions_path = write_dated_biogeography_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )

    assert "tree_is_time_scaled" in summary_path.read_text(encoding="utf-8")
    assert "age_before_present" in nodes_path.read_text(encoding="utf-8")
    assert "midpoint_age_before_present" in events_path.read_text(encoding="utf-8")
    assert "uncertainty_class" in bins_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")
