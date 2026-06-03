from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.regression import (
    summarize_phylogenetic_anova,
    write_phylogenetic_anova_exclusion_table,
    write_phylogenetic_anova_group_table,
    write_phylogenetic_anova_pairwise_table,
    write_phylogenetic_anova_simulation_table,
    write_phylogenetic_anova_summary_table,
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


def test_summarize_phylogenetic_anova_reports_seeded_group_effect() -> None:
    report = summarize_phylogenetic_anova(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_anova.tsv"),
        response="trait_value",
        group="habitat",
        simulations=32,
        seed=11,
    )
    replay = summarize_phylogenetic_anova(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_anova.tsv"),
        response="trait_value",
        group="habitat",
        simulations=32,
        seed=11,
    )

    assert report.tree_taxon_count == 6
    assert report.analyzed_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.analyzed_taxon_count == 6
    assert report.group_count == 2
    assert report.simulation_count == 32
    assert report.seed == 11
    assert report.pairwise_adjustment_method == "holm"
    assert report.f_statistic > 0.0
    assert 0.0 <= report.p_value <= 1.0
    assert math.isclose(report.f_statistic, replay.f_statistic)
    assert math.isclose(report.p_value, replay.p_value)
    assert [row.f_statistic for row in report.null_rows] == [
        row.f_statistic for row in replay.null_rows
    ]
    assert {row.group: row.taxon_count for row in report.group_rows} == {
        "desert": 2,
        "forest": 4,
    }
    assert report.low_sample_group_count == 1
    assert any("fewer than three taxa" in warning for warning in report.warnings)
    assert report.pairwise_rows[0].left_group == "desert"
    assert report.pairwise_rows[0].right_group == "forest"


def test_summarize_phylogenetic_anova_reports_missing_group_and_extra_taxon() -> None:
    report = summarize_phylogenetic_anova(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_anova_missing.tsv"),
        response="trait_value",
        group="habitat",
        simulations=16,
        seed=5,
    )

    assert report.analyzed_taxa == ["A", "B", "C", "D", "E"]
    assert {row.taxon: row.reason for row in report.excluded_taxa} == {
        "F": "missing_value",
        "G": "absent_from_tree",
    }
    assert {row.group: row.taxon_count for row in report.group_rows} == {
        "desert": 2,
        "forest": 3,
    }


def test_phylogenetic_anova_writers_emit_review_ledgers(tmp_path: Path) -> None:
    report = summarize_phylogenetic_anova(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_phylogenetic_anova.tsv"),
        response="trait_value",
        group="habitat",
        simulations=16,
        seed=7,
    )
    summary_out = tmp_path / "phylogenetic-anova-summary.tsv"
    groups_out = tmp_path / "phylogenetic-anova-groups.tsv"
    pairwise_out = tmp_path / "phylogenetic-anova-pairwise.tsv"
    simulation_out = tmp_path / "phylogenetic-anova-simulations.tsv"
    excluded_out = tmp_path / "phylogenetic-anova-excluded.tsv"

    write_phylogenetic_anova_summary_table(summary_out, report)
    write_phylogenetic_anova_group_table(groups_out, report)
    write_phylogenetic_anova_pairwise_table(pairwise_out, report)
    write_phylogenetic_anova_simulation_table(simulation_out, report)
    write_phylogenetic_anova_exclusion_table(excluded_out, report)

    assert (
        summary_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("response\tgroup\ttaxon_column")
    )
    assert (
        groups_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("group\ttaxon_count\ttaxa\tmean")
    )
    assert (
        pairwise_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("left_group\tright_group\tleft_taxon_count")
    )
    assert (
        simulation_out.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("simulation_index\tf_statistic\tat_or_above_observed")
    )
    assert excluded_out.read_text(encoding="utf-8").splitlines() == [
        "taxon\treason\tdetails"
    ]
