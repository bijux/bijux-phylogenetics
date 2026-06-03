from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.discrete.review import (
    summarize_ordered_discrete_reconstruction,
    summarize_ordered_discrete_report,
    write_ordered_discrete_fit_table,
    write_ordered_discrete_node_table,
    write_ordered_discrete_summary_table,
    write_ordered_discrete_transition_table,
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


def test_discrete_reconstruction_records_ordered_transition_restrictions() -> None:
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        state_ordering="ordered",
        ordered_states=["north", "south", "island"],
    )

    north_to_island = next(
        row
        for row in report.transition_rate_rows
        if row.source_state == "north" and row.target_state == "island"
    )
    north_to_south = next(
        row
        for row in report.transition_rate_rows
        if row.source_state == "north" and row.target_state == "south"
    )

    assert report.log_likelihood is not None
    assert report.parameter_count == 1
    assert report.aic is not None
    assert north_to_island.transition_allowed is False
    assert north_to_island.rate == 0.0
    assert north_to_south.transition_allowed is True
    assert north_to_south.rate > 0.0


def test_summarize_ordered_discrete_reconstruction_compares_ordered_and_unordered() -> (
    None
):
    report = summarize_ordered_discrete_reconstruction(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        ordered_states=["north", "south", "island"],
    )
    summary = summarize_ordered_discrete_report(report)

    assert [row.ordering_mode for row in report.fit_rows] == [
        "ordered",
        "unordered",
    ]
    assert summary.state_count == 3
    assert summary.ordered_parameter_count == 1
    assert summary.unordered_parameter_count == 1
    assert summary.ordered_log_likelihood != summary.unordered_log_likelihood
    assert summary.preferred_ordering == "ordered"
    assert summary.restricted_transition_count == 2
    assert len(report.node_rows) == 3
    north_to_island = next(
        row
        for row in report.transition_rows
        if row.source_state == "north" and row.target_state == "island"
    )
    assert north_to_island.ordered_transition_allowed is False
    assert north_to_island.unordered_transition_allowed is True
    assert north_to_island.ordered_rate == 0.0


def test_write_ordered_discrete_tables_emit_expected_ledgers(tmp_path: Path) -> None:
    report = summarize_ordered_discrete_reconstruction(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        ordered_states=["north", "south", "island"],
    )
    summary_path = tmp_path / "ordered-discrete-summary.tsv"
    fit_path = tmp_path / "ordered-discrete-fits.tsv"
    node_path = tmp_path / "ordered-discrete-nodes.tsv"
    transition_path = tmp_path / "ordered-discrete-transitions.tsv"

    write_ordered_discrete_summary_table(summary_path, report)
    write_ordered_discrete_fit_table(fit_path, report)
    write_ordered_discrete_node_table(node_path, report)
    write_ordered_discrete_transition_table(transition_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    fit_rows = fit_path.read_text(encoding="utf-8").splitlines()
    node_rows = node_path.read_text(encoding="utf-8").splitlines()
    transition_rows = transition_path.read_text(encoding="utf-8").splitlines()

    assert summary_rows[0].startswith(
        "trait\ttaxon_column\tmodel\tanalyzed_taxon_count\tstate_count"
    )
    assert fit_rows[0].startswith(
        "ordering_mode\tmodel\tstate_ordering\tordered_states"
    )
    assert node_rows[0].startswith(
        "node\tdescendant_taxa\tordered_state\tunordered_state"
    )
    assert transition_rows[0].startswith(
        "source_state\ttarget_state\tstep_distance\tordered_transition_allowed"
    )
    assert len(fit_rows) == 3
