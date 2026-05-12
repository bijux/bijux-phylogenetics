from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.biogeography import (
    summarize_geographic_state_model,
    write_geographic_exclusion_table,
    write_geographic_region_probability_table,
    write_geographic_state_summary_table,
    write_geographic_transition_event_table,
    write_geographic_transition_rate_table,
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


def test_summarize_geographic_state_model_supports_er_sym_and_ard_aliases() -> None:
    for model in ("er", "sym", "ard"):
        report = summarize_geographic_state_model(
            fixture("example_tree.nwk"),
            fixture("example_traits_geography.tsv"),
            trait="region",
            model=model,
        )

        assert report.model == model
        assert report.summary.model == model
        assert report.node_rows
        assert report.transition_rate_rows
        assert report.transition_event_rows


def test_summarize_geographic_state_model_tracks_explicit_exclusions(
    tmp_path: Path,
) -> None:
    traits_path = tmp_path / "geography.tsv"
    traits_path.write_text(
        "taxon\tregion\nA\tnorth\nB\tsouth\nC\tisland\nD\tsouth\nGhost\tisland\n",
        encoding="utf-8",
    )

    report = summarize_geographic_state_model(
        fixture("example_tree.nwk"),
        traits_path,
        trait="region",
        model="ard",
    )

    assert report.summary.excluded_taxon_count == 1
    assert report.exclusion_rows[0].taxon == "Ghost"
    assert report.exclusion_rows[0].reason == "taxon-not-in-tree"


def test_write_geographic_state_tables_emit_expected_ledgers(tmp_path: Path) -> None:
    report = summarize_geographic_state_model(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="ard",
    )

    summary_path = write_geographic_state_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    nodes_path = write_geographic_region_probability_table(
        tmp_path / "nodes.tsv",
        report,
    )
    rates_path = write_geographic_transition_rate_table(
        tmp_path / "rates.tsv",
        report,
    )
    events_path = write_geographic_transition_event_table(
        tmp_path / "events.tsv",
        report,
    )
    exclusions_path = write_geographic_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )

    assert "root_region" in summary_path.read_text(encoding="utf-8")
    assert "most_likely_region" in nodes_path.read_text(encoding="utf-8")
    assert "source_region" in rates_path.read_text(encoding="utf-8")
    assert "strongly_supported" in events_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")
