from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.tree_set import (
    summarize_continuous_ancestral_tree_set,
    summarize_continuous_ancestral_tree_set_report,
    summarize_discrete_ancestral_tree_set,
    summarize_discrete_ancestral_tree_set_report,
    write_ancestral_tree_set_exclusion_table,
    write_ancestral_tree_set_tree_table,
    write_continuous_ancestral_tree_set_clade_table,
    write_continuous_ancestral_tree_set_node_table,
    write_continuous_ancestral_tree_set_summary_table,
    write_discrete_ancestral_tree_set_clade_table,
    write_discrete_ancestral_tree_set_node_table,
    write_discrete_ancestral_tree_set_summary_table,
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


def test_summarize_continuous_ancestral_tree_set_reports_clade_distributions() -> None:
    report = summarize_continuous_ancestral_tree_set(
        fixture("example_posterior_tree_set_six_taxa.nwk"),
        fixture("example_traits_clade_summary.tsv"),
        trait="body_mass",
    )

    summary = summarize_continuous_ancestral_tree_set_report(report)
    root = next(
        row for row in report.clade_summaries if row.clade_id == "A|B|C|D|E|F"
    )

    assert report.total_tree_count == 5
    assert report.burnin_tree_count == 0
    assert report.kept_tree_count == 5
    assert report.shared_tree_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.analysis_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.rooted_topology_count == 5
    assert report.unrooted_topology_count == 4
    assert len(report.tree_rows) == 5
    assert len(report.node_rows) == 25
    assert len(report.clade_summaries) == 14
    assert root.tree_presence_count == 5
    assert root.tree_presence_fraction == 1.0
    assert root.mean_estimate > 5.0
    assert root.stability_class == "within_tree_uncertainty"
    assert summary.clade_summary_count == 14
    assert summary.unstable_clade_count == 14
    assert summary.top_unstable_clade is not None
    assert (
        "one or more comparable ancestral clades are absent from some retained trees"
        in report.warnings
    )


def test_summarize_discrete_ancestral_tree_set_reports_state_distributions() -> None:
    report = summarize_discrete_ancestral_tree_set(
        fixture("example_posterior_tree_set_six_taxa.nwk"),
        fixture("example_traits_clade_summary.tsv"),
        trait="habitat",
        model="equal-rates",
    )

    summary = summarize_discrete_ancestral_tree_set_report(report)
    island_clade = next(row for row in report.clade_summaries if row.clade_id == "E|F")
    root = next(
        row for row in report.clade_summaries if row.clade_id == "A|B|C|D|E|F"
    )

    assert report.total_tree_count == 5
    assert report.kept_tree_count == 5
    assert report.ordered_states == []
    assert report.shared_tree_taxa == ["A", "B", "C", "D", "E", "F"]
    assert len(report.tree_rows) == 5
    assert len(report.node_rows) == 25
    assert len(report.clade_summaries) == 14
    assert island_clade.tree_presence_count == 3
    assert island_clade.dominant_state == "island"
    assert island_clade.state_distribution == {"island": 3}
    assert island_clade.stability_class == "topology_sensitive"
    assert root.dominant_state == "forest"
    assert root.dominant_state_fraction == 1.0
    assert root.stability_class == "low_confidence"
    assert summary.clade_summary_count == 14
    assert summary.unstable_clade_count == 14
    assert (
        "one or more discrete ancestral clades change state or support profile across retained trees"
        in report.warnings
    )


def test_summarize_continuous_ancestral_tree_set_applies_burnin_and_exclusions(
    tmp_path: Path,
) -> None:
    table_path = tmp_path / "traits.tsv"
    table_path.write_text(
        "\n".join(
            [
                "taxon\tbody_mass\thabitat",
                "A\t1.0\tforest",
                "B\t2.0\tforest",
                "C\t2.5\tforest",
                "D\t3.0\ttundra",
                "E\t10.0\tisland",
                "F\t\tisland",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    report = summarize_continuous_ancestral_tree_set(
        fixture("example_posterior_tree_set_six_taxa.nwk"),
        table_path,
        trait="body_mass",
        burnin_fraction=0.2,
    )

    assert report.burnin_tree_count == 1
    assert report.kept_tree_count == 4
    assert sorted(report.analysis_taxa) == ["A", "B", "C", "D", "E"]
    assert [(row.taxon, row.reason) for row in report.exclusions] == [
        ("F", "missing_trait_value"),
    ]


def test_write_ancestral_tree_set_tables_write_expected_rows(tmp_path: Path) -> None:
    continuous_report = summarize_continuous_ancestral_tree_set(
        fixture("example_posterior_tree_set_six_taxa.nwk"),
        fixture("example_traits_clade_summary.tsv"),
        trait="body_mass",
    )
    discrete_report = summarize_discrete_ancestral_tree_set(
        fixture("example_posterior_tree_set_six_taxa.nwk"),
        fixture("example_traits_clade_summary.tsv"),
        trait="habitat",
        model="equal-rates",
    )

    continuous_summary_path = tmp_path / "continuous-summary.tsv"
    continuous_node_path = tmp_path / "continuous-nodes.tsv"
    continuous_clade_path = tmp_path / "continuous-clades.tsv"
    discrete_summary_path = tmp_path / "discrete-summary.tsv"
    discrete_node_path = tmp_path / "discrete-nodes.tsv"
    discrete_clade_path = tmp_path / "discrete-clades.tsv"
    tree_path = tmp_path / "trees.tsv"
    exclusion_path = tmp_path / "excluded.tsv"

    write_ancestral_tree_set_tree_table(tree_path, continuous_report)
    write_ancestral_tree_set_exclusion_table(exclusion_path, continuous_report)
    write_continuous_ancestral_tree_set_summary_table(
        continuous_summary_path,
        continuous_report,
    )
    write_continuous_ancestral_tree_set_node_table(
        continuous_node_path,
        continuous_report,
    )
    write_continuous_ancestral_tree_set_clade_table(
        continuous_clade_path,
        continuous_report,
    )
    write_discrete_ancestral_tree_set_summary_table(
        discrete_summary_path,
        discrete_report,
    )
    write_discrete_ancestral_tree_set_node_table(
        discrete_node_path,
        discrete_report,
    )
    write_discrete_ancestral_tree_set_clade_table(
        discrete_clade_path,
        discrete_report,
    )

    tree_rows = tree_path.read_text(encoding="utf-8").splitlines()
    exclusion_rows = exclusion_path.read_text(encoding="utf-8").splitlines()
    continuous_summary_rows = continuous_summary_path.read_text(
        encoding="utf-8"
    ).splitlines()
    continuous_node_rows = continuous_node_path.read_text(encoding="utf-8").splitlines()
    continuous_clade_rows = continuous_clade_path.read_text(
        encoding="utf-8"
    ).splitlines()
    discrete_summary_rows = discrete_summary_path.read_text(
        encoding="utf-8"
    ).splitlines()
    discrete_node_rows = discrete_node_path.read_text(encoding="utf-8").splitlines()
    discrete_clade_rows = discrete_clade_path.read_text(encoding="utf-8").splitlines()

    assert tree_rows[0].startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id"
    )
    assert exclusion_rows == ["taxon\treason"]
    assert continuous_summary_rows[0].startswith("trait\ttaxon_column\tmodel\talpha")
    assert continuous_node_rows[0].startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id\tunrooted_topology_id\tclade_id"
    )
    assert continuous_clade_rows[0].startswith(
        "clade_id\tclade_taxa\ttree_presence_count\ttree_presence_fraction"
    )
    assert discrete_summary_rows[0].startswith(
        "trait\ttaxon_column\tmodel\tstate_ordering"
    )
    assert discrete_node_rows[0].startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id\tunrooted_topology_id\tclade_id"
    )
    assert discrete_clade_rows[0].startswith(
        "clade_id\tclade_taxa\ttree_presence_count\ttree_presence_fraction"
    )
    assert len(tree_rows) == 6
    assert len(continuous_summary_rows) == 2
    assert len(continuous_node_rows) == 26
    assert len(continuous_clade_rows) == 15
    assert len(discrete_summary_rows) == 2
    assert len(discrete_node_rows) == 26
    assert len(discrete_clade_rows) == 15
