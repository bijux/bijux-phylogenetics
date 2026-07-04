from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.biogeography.migration import (
    summarize_geographic_migration_event_tree_set,
    summarize_geographic_migration_events,
    write_geographic_migration_event_summary_table,
    write_geographic_migration_event_table,
    write_geographic_migration_exclusion_table,
    write_geographic_migration_tree_set_event_summary_table,
    write_geographic_migration_tree_set_event_table,
    write_geographic_migration_tree_set_exclusion_table,
    write_geographic_migration_tree_set_summary_table,
    write_geographic_migration_tree_set_tree_table,
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


def test_summarize_geographic_migration_events_reports_branch_depths() -> None:
    report = summarize_geographic_migration_events(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="ard",
    )

    assert report.summary.event_count == 2
    assert report.summary.tree_depth == 0.3
    assert [row.source_region for row in report.event_rows] == ["north", "south"]
    assert [row.target_region for row in report.event_rows] == ["south", "island"]
    assert report.event_rows[0].branch_id == "C|D"
    assert report.event_rows[0].midpoint_depth == 0.05
    assert report.event_rows[1].branch_id == "D"
    assert report.event_rows[1].midpoint_depth == 0.2


def test_summarize_geographic_migration_event_tree_set_reports_topology_sensitive_events() -> (
    None
):
    report = summarize_geographic_migration_event_tree_set(
        fixture("example_tree_set_left.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="ard",
    )

    assert report.summary.kept_tree_count == 3
    assert report.summary.event_row_count == len(report.event_rows)
    assert report.summary.event_summary_count == len(report.event_summaries)
    assert report.summary.topology_sensitive_event_count >= 1
    assert any(
        row.stability_class == "topology_sensitive" for row in report.event_summaries
    )
    assert any("topology-sensitive" in warning for warning in report.warnings)


def test_summarize_geographic_migration_event_tree_set_supports_burnin() -> None:
    report = summarize_geographic_migration_event_tree_set(
        fixture("example_tree_set_left.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="sym",
        burnin_fraction=1 / 3,
    )

    assert report.summary.total_tree_count == 3
    assert report.summary.burnin_tree_count == 1
    assert report.summary.kept_tree_count == 2
    assert all(row.post_burnin_index >= 1 for row in report.tree_rows)


def test_write_geographic_migration_event_tables_emit_expected_ledgers(
    tmp_path: Path,
) -> None:
    report = summarize_geographic_migration_events(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="er",
    )

    summary_path = write_geographic_migration_event_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    events_path = write_geographic_migration_event_table(
        tmp_path / "events.tsv",
        report,
    )
    exclusions_path = write_geographic_migration_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )

    assert "mean_event_support" in summary_path.read_text(encoding="utf-8")
    assert "midpoint_depth" in events_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")


def test_write_geographic_migration_tree_set_tables_emit_expected_ledgers(
    tmp_path: Path,
) -> None:
    report = summarize_geographic_migration_event_tree_set(
        fixture("example_tree_set_left.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="er",
    )

    summary_path = write_geographic_migration_tree_set_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    trees_path = write_geographic_migration_tree_set_tree_table(
        tmp_path / "trees.tsv",
        report,
    )
    events_path = write_geographic_migration_tree_set_event_table(
        tmp_path / "events.tsv",
        report,
    )
    event_summaries_path = write_geographic_migration_tree_set_event_summary_table(
        tmp_path / "event-summaries.tsv",
        report,
    )
    exclusions_path = write_geographic_migration_tree_set_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )

    assert "event_summary_count" in summary_path.read_text(encoding="utf-8")
    assert "rooted_topology_id" in trees_path.read_text(encoding="utf-8")
    assert "midpoint_depth" in events_path.read_text(encoding="utf-8")
    assert "tree_presence_fraction" in event_summaries_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")
