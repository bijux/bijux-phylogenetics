from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.clades.stability import (
    analyze_comparative_clade_stability,
    write_comparative_clade_coefficient_change_table,
    write_comparative_clade_stability_table,
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


def test_analyze_comparative_clade_stability_ranks_influential_pgls_clades() -> None:
    report = analyze_comparative_clade_stability(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        formula="response_growth ~ predictor_two",
        lambda_value=0.0,
    )

    assert report.model_family == "pgls"
    assert report.baseline_term_count == 2
    assert report.minimum_major_clade_size == 2
    assert report.candidate_clade_count == 4
    assert report.blocked_clade_count == 1
    assert report.influential_clades == ["A|B", "C|D", "E|F"]

    ranked = sorted(report.clade_rows, key=lambda row: (row.rank or 999, row.clade_id))
    assert ranked[0].clade_id == "A|B"
    assert ranked[0].sign_changed_term_count == 1
    assert ranked[0].coefficient_comparison_count == 2
    assert ranked[0].missing_baseline_term_count == 0

    blocked = next(row for row in report.clade_rows if row.clade_id == "A|B|C|D")
    assert blocked.fit_status == "blocked"
    assert blocked.rank == 0
    assert blocked.missing_baseline_term_count == 2
    assert "insufficient for comparative refitting" in (blocked.blocker or "")

    predictor_two = next(
        row
        for row in report.coefficient_rows
        if row.clade_id == "A|B" and row.term == "predictor_two"
    )
    assert predictor_two.sign_changed is True
    assert predictor_two.significance_changed is False
    assert predictor_two.delta_estimate < 0.0


def test_analyze_comparative_clade_stability_auto_detects_logistic_family() -> None:
    report = analyze_comparative_clade_stability(
        fixture("example_tree_eight_taxa.nwk"),
        fixture("example_traits_phylogenetic_logistic_model_selection.tsv"),
        formula="presence ~ body_size",
        lambda_value=1.0,
    )

    assert report.model_family == "logistic"
    assert report.candidate_clade_count == 6
    assert report.blocked_clade_count == 2
    assert report.baseline_term_count == 2

    blocked_ids = {
        row.clade_id for row in report.clade_rows if row.fit_status == "blocked"
    }
    assert blocked_ids == {"A|B|C|D", "E|F|G|H"}
    for row in report.clade_rows:
        if row.fit_status == "blocked":
            assert "requires at least one success and one failure" in (
                row.blocker or ""
            )


def test_write_comparative_clade_stability_tables_write_expected_rows(
    tmp_path: Path,
) -> None:
    report = analyze_comparative_clade_stability(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        formula="response_growth ~ predictor_two",
        lambda_value=0.0,
    )
    summary_out = tmp_path / "comparative-clade-stability.tsv"
    coefficients_out = tmp_path / "comparative-clade-coefficients.tsv"

    write_comparative_clade_stability_table(summary_out, report)
    write_comparative_clade_coefficient_change_table(coefficients_out, report)

    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    coefficient_rows = coefficients_out.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("clade_id\tnode_label\tdropped_taxon_count")
    assert coefficient_rows[0].startswith(
        "clade_id\tnode_label\tterm\tbaseline_estimate"
    )
    assert len(summary_rows) == 5
    assert len(coefficient_rows) == 7
