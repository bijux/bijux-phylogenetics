from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.clades.traits import (
    summarize_clade_traits,
    write_clade_trait_clade_table,
    write_clade_trait_exclusion_table,
    write_clade_trait_summary_table,
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


def test_summarize_clade_traits_reports_continuous_clade_statistics() -> None:
    report = summarize_clade_traits(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_clade_summary.tsv"),
        trait="body_mass",
    )

    assert report.trait_kind == "continuous"
    assert report.tree_taxon_count == 6
    assert report.analyzed_taxa == ["A", "B", "C", "D", "E", "F"]
    assert len(report.clade_rows) == 4
    assert report.top_exceptional_clade == "E|F"
    assert report.exceptional_clades == ["E|F"]
    assert math.isclose(report.baseline_mean or 0.0, 4.916666666666667)
    assert math.isclose(report.baseline_median or 0.0, 2.75)
    assert math.isclose(report.baseline_range_width or 0.0, 10.0)

    top_row = report.clade_rows[0]
    assert top_row.clade_id == "E|F"
    assert math.isclose(top_row.mean or 0.0, 10.5)
    assert math.isclose(top_row.median or 0.0, 10.5)
    assert math.isclose(top_row.minimum or 0.0, 10.0)
    assert math.isclose(top_row.maximum or 0.0, 11.0)
    assert math.isclose(top_row.range_width or 0.0, 1.0)
    assert top_row.exceptional is True
    assert top_row.rank == 1


def test_summarize_clade_traits_reports_categorical_clade_distributions() -> None:
    report = summarize_clade_traits(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_clade_summary.tsv"),
        trait="habitat",
    )

    assert report.trait_kind == "categorical"
    assert report.top_exceptional_clade == "E|F"
    assert report.baseline_dominant_state == "forest"
    assert math.isclose(report.baseline_dominant_state_fraction or 0.0, 0.5)

    top_row = report.clade_rows[0]
    assert top_row.clade_id == "E|F"
    assert top_row.dominant_state == "island"
    assert top_row.dominant_state_count == 2
    assert math.isclose(top_row.dominant_state_fraction or 0.0, 1.0)
    assert top_row.distinct_state_count == 1
    assert [f"{row.state}={row.count}" for row in top_row.state_counts] == ["island=2"]
    assert math.isclose(top_row.distribution_shift or 0.0, 2.0 / 3.0)
    assert top_row.exceptional is True
    assert top_row.rank == 1


def test_summarize_clade_traits_tracks_excluded_taxa() -> None:
    report = summarize_clade_traits(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_continuous_evolution_missing.tsv"),
        trait="response_growth",
        trait_kind="continuous",
    )

    assert report.analyzed_taxa == ["A", "D", "E", "F"]
    assert {row.taxon: row.reason for row in report.excluded_taxa} == {
        "B": "missing_trait_value",
        "C": "non_numeric_trait_value",
        "G": "absent_from_tree",
    }


def test_clade_trait_writers_emit_review_ledgers(tmp_path: Path) -> None:
    report = summarize_clade_traits(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_clade_summary.tsv"),
        trait="habitat",
    )
    summary_out = tmp_path / "clade-traits-summary.tsv"
    clades_out = tmp_path / "clade-traits.tsv"
    excluded_out = tmp_path / "clade-traits-excluded.tsv"

    write_clade_trait_summary_table(summary_out, report)
    write_clade_trait_clade_table(clades_out, report)
    write_clade_trait_exclusion_table(excluded_out, report)

    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    clade_rows = clades_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("trait\ttaxon_column\ttrait_kind")
    assert clade_rows[0].startswith("clade_id\tnode_label\ttrait_kind")
    assert "E|F" in clade_rows[1]
    assert excluded_rows == ["taxon\treason"]
