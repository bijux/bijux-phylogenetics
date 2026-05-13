from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.biogeography import (
    summarize_geographic_sampling_bias,
    write_geographic_sampling_bias_exclusion_table,
    write_geographic_sampling_bias_node_table,
    write_geographic_sampling_bias_summary_table,
    write_geographic_sampling_bias_transition_table,
    write_geographic_sampling_count_table,
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


def test_summarize_geographic_sampling_bias_reports_region_counts_and_flags_dominance() -> (
    None
):
    report = summarize_geographic_sampling_bias(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_geography_dominated.tsv"),
        trait="region",
        model="ard",
    )

    assert report.summary.model == "ard"
    assert report.summary.weighting_mode == "inverse-frequency"
    assert report.summary.region_dominated is True
    assert report.summary.dominant_region == "north"
    assert report.summary.dominant_region_fraction == 0.833333333333333
    assert report.summary.weighted_dominant_region_fraction == 0.5
    assert len(report.count_rows) == 2
    north_row = next(row for row in report.count_rows if row.region == "north")
    assert north_row.sample_count == 5
    assert north_row.applied_weight == 0.6
    assert north_row.weighted_sample_fraction == 0.5


def test_summarize_geographic_sampling_bias_supports_explicit_region_weights() -> None:
    report = summarize_geographic_sampling_bias(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_geography_biased.tsv"),
        trait="region",
        model="ard",
        weights_path=fixture("example_geographic_region_weights.tsv"),
    )

    assert report.summary.weighting_mode == "explicit"
    assert report.summary.root_region_unweighted == "north"
    assert report.summary.root_region_weighted == "south"
    assert report.summary.root_region_changed is True
    assert report.summary.changed_internal_node_count >= 1
    assert any(row.changed for row in report.node_rows if row.is_root)
    assert any(row.changed_by_weighting for row in report.transition_rows)


def test_write_geographic_sampling_bias_tables_emit_expected_ledgers(
    tmp_path: Path,
) -> None:
    report = summarize_geographic_sampling_bias(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_geography_biased.tsv"),
        trait="region",
        model="ard",
        weights_path=fixture("example_geographic_region_weights.tsv"),
    )

    summary_path = write_geographic_sampling_bias_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    counts_path = write_geographic_sampling_count_table(
        tmp_path / "counts.tsv",
        report,
    )
    nodes_path = write_geographic_sampling_bias_node_table(
        tmp_path / "nodes.tsv",
        report,
    )
    transitions_path = write_geographic_sampling_bias_transition_table(
        tmp_path / "transitions.tsv",
        report,
    )
    exclusions_path = write_geographic_sampling_bias_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )

    assert "weighting_mode" in summary_path.read_text(encoding="utf-8")
    assert "weighted_sample_fraction" in counts_path.read_text(encoding="utf-8")
    assert "weighted_region_probabilities" in nodes_path.read_text(encoding="utf-8")
    assert "changed_by_weighting" in transitions_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")
