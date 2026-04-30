from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.reporting import (
    build_comparative_method_report,
    build_trait_influence_report,
    compare_comparative_results_across_pruning,
    compare_comparative_results_across_trees,
    write_comparative_method_report,
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
    assert "causal interpretation is not warranted from comparative association alone" in report.snapshot.limitations


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
    assert "maturity" in out_path.read_text(encoding="utf-8")
