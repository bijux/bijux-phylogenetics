from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.reporting import (
    build_comparative_method_report,
    build_comparative_methods_summary_text,
    build_trait_influence_report,
    compare_comparative_results_across_pruning,
    compare_comparative_results_across_trees,
    write_comparative_method_report,
    write_comparative_methods_summary_text,
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


def test_build_comparative_method_report_returns_audit_rows_and_limitations() -> None:
    report = build_comparative_method_report(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=0.0,
    )
    assert len(report.snapshot.audit_rows) == 3
    assert report.snapshot.maturity.reference_validation_passed is True
    assert report.snapshot.pgls_inputs.formula_audit.encoded_columns
    assert (
        "causal interpretation is not warranted from comparative association alone"
        in report.snapshot.limitations
    )


def test_build_trait_influence_report_combines_predictor_and_taxon_rankings() -> None:
    report = build_trait_influence_report(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=0.0,
    )
    assert len(report.predictor_rows) == 1
    assert len(report.taxon_rows) == 4
    assert report.top_predictor_terms == ["predictor_one"]
    assert len(report.top_taxa) == 3
    assert report.taxon_rows[0].taxon in {"A", "D", "B"}


def test_compare_comparative_results_across_trees_reports_metric_deltas() -> None:
    report = compare_comparative_results_across_trees(
        fixture("example_tree.nwk"),
        fixture("example_tree_topology_diff.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=0.0,
    )
    assert report.coefficient_deltas
    assert report.left_selected_model in {"brownian", "ou"}
    assert report.right_selected_model in {"brownian", "ou"}
    assert isinstance(report.conclusion_changed, bool)
    assert report.sign_changed_terms == []


def test_compare_comparative_results_across_pruning_tracks_dropped_taxa() -> None:
    report = compare_comparative_results_across_pruning(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_interaction.tsv"),
        formula="response ~ predictor_one + habitat",
        drop_taxa=["F"],
        lambda_value=0.0,
    )
    assert report.dropped_taxa == ["F"]
    assert len(report.pruned_taxa) == 5
    assert isinstance(report.conclusion_changed, bool)


def test_write_comparative_method_report_writes_html(tmp_path: Path) -> None:
    report = build_comparative_method_report(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=0.0,
    )
    out_path = tmp_path / "comparative-report.html"
    write_comparative_method_report(out_path, report)
    assert out_path.exists()
    assert "Bijux Comparative Method Report" in out_path.read_text(encoding="utf-8")
    assert "formula-audit" in out_path.read_text(encoding="utf-8")
    assert "maturity" in out_path.read_text(encoding="utf-8")


def test_build_comparative_methods_summary_text_reports_pruning_model_and_diagnostics() -> (
    None
):
    report = build_comparative_method_report(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_interaction.tsv"),
        formula="response ~ predictor_one + habitat",
        lambda_value=0.0,
    )

    text = build_comparative_methods_summary_text(report)

    assert "Comparative Analysis Methods Summary" in text
    assert "- tree taxon count: `6`" in text
    assert "- analysis taxa retained after overlap and numeric pruning: `6`" in text
    assert "- formula-level excluded taxa: `0`" in text
    assert "- categorical predictors: `habitat`" in text
    assert "- predictor terms: `predictor_one`, `habitat`" in text
    assert "- selected comparative process model:" in text
    assert "- PGLS R-squared:" in text
    assert "- diagnostic and reviewer warnings:" in text


def test_write_comparative_methods_summary_text_writes_markdown(tmp_path: Path) -> None:
    report = build_comparative_method_report(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=0.0,
    )
    output_path = tmp_path / "comparative-methods-summary.md"

    result = write_comparative_methods_summary_text(output_path, report)

    assert result.output_path == output_path
    assert result.selected_model in {"brownian", "ou"}
    assert result.predictor_count == 1
    assert result.analysis_taxa == 4
    assert result.excluded_taxa == 0
    assert "Comparative Analysis Methods Summary" in result.text
    assert output_path.read_text(encoding="utf-8") == result.text
