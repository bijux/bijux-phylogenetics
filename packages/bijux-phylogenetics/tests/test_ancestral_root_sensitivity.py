from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.sensitivity import (
    summarize_ancestral_root_sensitivity,
    summarize_ancestral_root_sensitivity_report,
    write_ancestral_root_assumption_table,
    write_ancestral_root_sensitivity_node_table,
    write_ancestral_root_sensitivity_summary_table,
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


def test_discrete_reconstruction_supports_empirical_and_fixed_root_priors() -> None:
    equal_report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        root_prior_mode="equal",
    )
    empirical_report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        root_prior_mode="empirical",
    )
    fixed_report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        root_prior_mode="fixed",
        fixed_root_state="island",
    )

    equal_root = next(
        estimate for estimate in equal_report.estimates if estimate.node == "A|B|C|D"
    )
    empirical_root = next(
        estimate
        for estimate in empirical_report.estimates
        if estimate.node == "A|B|C|D"
    )
    fixed_root = next(
        estimate for estimate in fixed_report.estimates if estimate.node == "A|B|C|D"
    )

    assert equal_root.most_likely_state == "north"
    assert empirical_root.most_likely_state == "north"
    assert empirical_root.confidence > equal_root.confidence
    assert fixed_root.most_likely_state == "island"
    assert fixed_root.confidence == 1.0


def test_summarize_ancestral_root_sensitivity_reports_state_and_support_changes() -> (
    None
):
    report = summarize_ancestral_root_sensitivity(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        fixed_root_state="island",
    )
    summary = summarize_ancestral_root_sensitivity_report(report)

    assert [row.assumption_id for row in report.assumption_rows] == [
        "equal_root_prior",
        "empirical_root_prior",
        "fixed_root_state",
    ]
    assert summary.assumption_count == 3
    assert summary.compared_node_count == 3
    assert summary.state_changed_node_count == 2
    assert summary.support_changed_node_count == 1
    assert summary.top_sensitive_node == "A|B|C|D"
    assert report.node_rows[0].stability_class == "root_sensitive_state"
    assert report.node_rows[1].stability_class == "root_sensitive_state"
    assert "equal_root_prior" in report.node_rows[0].assumption_states


def test_write_ancestral_root_sensitivity_tables_emit_expected_rows(
    tmp_path: Path,
) -> None:
    report = summarize_ancestral_root_sensitivity(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
        fixed_root_state="island",
    )
    summary_path = tmp_path / "ancestral-root-sensitivity-summary.tsv"
    assumptions_path = tmp_path / "ancestral-root-assumptions.tsv"
    nodes_path = tmp_path / "ancestral-root-nodes.tsv"

    write_ancestral_root_sensitivity_summary_table(summary_path, report)
    write_ancestral_root_assumption_table(assumptions_path, report)
    write_ancestral_root_sensitivity_node_table(nodes_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    assumption_rows = assumptions_path.read_text(encoding="utf-8").splitlines()
    node_rows = nodes_path.read_text(encoding="utf-8").splitlines()

    assert summary_rows[0].startswith(
        "trait\ttaxon_column\tmodel\tstate_ordering\tanalyzed_taxon_count"
    )
    assert assumption_rows[0].startswith(
        "assumption_id\troot_prior_mode\tfixed_root_state\troot_prior_distribution"
    )
    assert node_rows[0].startswith(
        "node\tdescendant_taxa\tassumption_states\tassumption_confidences"
    )
    assert len(assumption_rows) == 4
    assert len(node_rows) == 4
