from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.continuous import (
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.presentation.confidence_review import (
    build_continuous_ancestral_confidence_rows,
    build_continuous_ancestral_tree_set_confidence_rows,
    build_discrete_ancestral_confidence_rows,
    build_discrete_ancestral_tree_set_confidence_rows,
    summarize_continuous_ancestral_confidence,
    summarize_continuous_ancestral_tree_set_confidence,
    summarize_discrete_ancestral_confidence,
    summarize_discrete_ancestral_tree_set_confidence,
    write_ancestral_confidence_summary_table,
    write_continuous_ancestral_confidence_table,
    write_discrete_ancestral_confidence_table,
)
from bijux_phylogenetics.ancestral.tree_set import (
    summarize_continuous_ancestral_tree_set,
    summarize_discrete_ancestral_tree_set,
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


def test_build_continuous_ancestral_confidence_rows_ranks_internal_nodes() -> None:
    report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )

    rows = build_continuous_ancestral_confidence_rows(report)
    summary = summarize_continuous_ancestral_confidence(report)

    assert [row.uncertainty_rank for row in rows] == [1, 2, 3]
    assert rows[0].node == "A|B|C|D"
    assert rows[0].uncertainty_score >= rows[1].uncertainty_score
    assert rows[-1].node == "A|B"
    assert all(row.confidence_class == "low" for row in rows)
    assert summary.top_uncertain_id == "A|B|C|D"
    assert summary.confidence_row_count == 3
    assert summary.low_confidence_count == 3


def test_build_discrete_ancestral_confidence_rows_reports_entropy_and_posterior() -> (
    None
):
    report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )

    rows = build_discrete_ancestral_confidence_rows(report)
    summary = summarize_discrete_ancestral_confidence(report)

    assert [row.uncertainty_rank for row in rows] == [1, 2, 3]
    assert rows[0].node == "C|D"
    assert rows[0].max_posterior_probability < 0.5
    assert rows[0].normalized_entropy > 0.9
    assert rows[1].node == "A|B|C|D"
    assert rows[2].confidence_class == "moderate"
    assert summary.top_uncertain_id == "C|D"
    assert summary.high_entropy_count == 2
    assert summary.low_confidence_count == 2


def test_build_ancestral_tree_set_confidence_rows_supports_continuous_and_discrete() -> (
    None
):
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

    continuous_rows = build_continuous_ancestral_tree_set_confidence_rows(
        continuous_report
    )
    discrete_rows = build_discrete_ancestral_tree_set_confidence_rows(discrete_report)
    continuous_summary = summarize_continuous_ancestral_tree_set_confidence(
        continuous_report
    )
    discrete_summary = summarize_discrete_ancestral_tree_set_confidence(discrete_report)

    assert continuous_rows[0].clade_id == "A|B|C|D|E|F"
    assert continuous_rows[0].uncertainty_rank == 1
    assert continuous_rows[0].confidence_class == "low"
    assert continuous_summary.kept_tree_count == 5
    assert continuous_summary.top_uncertain_id == "A|B|C|D|E|F"

    assert discrete_rows[0].tree_presence_fraction == 0.2
    assert discrete_rows[0].uncertainty_score >= discrete_rows[1].uncertainty_score
    assert discrete_rows[0].confidence_class == "low"
    assert discrete_summary.kept_tree_count == 5
    assert discrete_summary.low_confidence_count >= 10
    assert discrete_summary.top_uncertain_id is not None


def test_write_ancestral_confidence_tables_emit_ranked_ledgers(tmp_path: Path) -> None:
    continuous_report = reconstruct_continuous_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        trait="response",
    )
    discrete_report = reconstruct_discrete_ancestral_states(
        fixture("example_tree.nwk"),
        fixture("example_traits_geography.tsv"),
        trait="region",
        model="equal-rates",
    )
    summary_path = tmp_path / "ancestral-confidence-summary.tsv"
    continuous_path = tmp_path / "continuous-confidence.tsv"
    discrete_path = tmp_path / "discrete-confidence.tsv"

    write_ancestral_confidence_summary_table(
        summary_path,
        summarize_discrete_ancestral_confidence(discrete_report),
    )
    write_continuous_ancestral_confidence_table(continuous_path, continuous_report)
    write_discrete_ancestral_confidence_table(discrete_path, discrete_report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    continuous_rows = continuous_path.read_text(encoding="utf-8").splitlines()
    discrete_rows = discrete_path.read_text(encoding="utf-8").splitlines()

    assert summary_rows[0].startswith(
        "trait\ttaxon_column\tsource_kind\treconstruction_kind\ttarget_kind"
    )
    assert continuous_rows[0].startswith(
        "node\tnode_name\tdescendant_taxa\testimate\tstandard_error"
    )
    assert discrete_rows[0].startswith(
        "node\tnode_name\tdescendant_taxa\tmost_likely_state\tstate_set"
    )
    assert "\t1\tlow\ttrue" in continuous_rows[1]
    assert "normalized_entropy" in discrete_rows[0]
