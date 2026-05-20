from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.ecology import (
    summarize_niche_transitions,
    write_niche_state_node_table,
    write_niche_transition_branch_table,
    write_niche_transition_clade_table,
    write_niche_transition_count_table,
    write_niche_transition_exclusion_table,
    write_niche_transition_rate_table,
    write_niche_transition_summary_table,
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


@pytest.mark.slow
def test_summarize_niche_transitions_supports_er_sym_and_ard_aliases() -> None:
    for model in ("er", "sym", "ard"):
        report = summarize_niche_transitions(
            fixture("example_tree_six_taxa.nwk"),
            fixture("example_traits_ecological_niche.tsv"),
            trait="niche",
            model=model,
        )

        assert report.model == model
        assert report.summary.model == model
        assert report.node_rows
        assert report.rate_rows
        assert report.branch_rows
        assert report.count_rows
        assert report.clade_rows


def test_summarize_niche_transitions_reports_clade_specific_shift_burden() -> None:
    report = summarize_niche_transitions(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_ecological_niche.tsv"),
        trait="niche",
        model="er",
    )

    top_clade = report.clade_rows[0]

    assert report.summary.changed_branch_count >= 2
    assert report.summary.clade_shift_row_count == len(report.clade_rows)
    assert top_clade.rank == 1
    assert top_clade.changed_branch_count >= 1
    assert top_clade.shift_burden_score > 0.0


def test_summarize_niche_transitions_tracks_explicit_exclusions(
    tmp_path: Path,
) -> None:
    traits_path = tmp_path / "ecology.tsv"
    traits_path.write_text(
        (
            "taxon\tniche\n"
            "A\tforest\n"
            "B\tgrassland\n"
            "C\tdesert\n"
            "D\tdesert\n"
            "E\tmarine\n"
            "F\tmarine\n"
            "Ghost\tmarine\n"
            "Dropped\t\n"
        ),
        encoding="utf-8",
    )

    report = summarize_niche_transitions(
        fixture("example_tree_six_taxa.nwk"),
        traits_path,
        trait="niche",
        model="er",
    )

    assert report.summary.excluded_taxon_count == 2
    assert {(row.taxon, row.reason) for row in report.exclusion_rows} == {
        ("Ghost", "taxon-not-in-tree"),
        ("Dropped", "missing-state"),
    }


def test_write_niche_transition_tables_emit_expected_ledgers(tmp_path: Path) -> None:
    report = summarize_niche_transitions(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_ecological_niche.tsv"),
        trait="niche",
        model="er",
    )

    summary_path = write_niche_transition_summary_table(
        tmp_path / "summary.tsv", report
    )
    nodes_path = write_niche_state_node_table(tmp_path / "nodes.tsv", report)
    rates_path = write_niche_transition_rate_table(tmp_path / "rates.tsv", report)
    branches_path = write_niche_transition_branch_table(
        tmp_path / "branches.tsv", report
    )
    counts_path = write_niche_transition_count_table(tmp_path / "counts.tsv", report)
    clades_path = write_niche_transition_clade_table(tmp_path / "clades.tsv", report)
    exclusions_path = write_niche_transition_exclusion_table(
        tmp_path / "exclusions.tsv",
        report,
    )

    assert "root_niche" in summary_path.read_text(encoding="utf-8")
    assert "most_likely_niche" in nodes_path.read_text(encoding="utf-8")
    assert "source_niche" in rates_path.read_text(encoding="utf-8")
    assert "certainty_class" in branches_path.read_text(encoding="utf-8")
    assert "certain_transition_count" in counts_path.read_text(encoding="utf-8")
    assert "shift_burden_score" in clades_path.read_text(encoding="utf-8")
    assert "reason" in exclusions_path.read_text(encoding="utf-8")
