from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.pgls.posterior_tree import (
    run_posterior_tree_pgls,
    write_posterior_tree_pgls_coefficient_table,
    write_posterior_tree_pgls_summary_table,
    write_posterior_tree_pgls_tree_table,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected", "mrbayes")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_run_posterior_tree_pgls_summarizes_coefficient_stability() -> None:
    report = run_posterior_tree_pgls(
        fixture("example_posterior_tree_set_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        formula="response_growth ~ predictor_two",
        lambda_value="estimate",
        significance_threshold=0.1,
    )

    assert report.response == "response_growth"
    assert report.lambda_mode == "estimate"
    assert report.total_tree_count == 5
    assert report.burnin_tree_count == 0
    assert report.kept_tree_count == 5
    assert report.rooted_topology_count == 5
    assert report.unrooted_topology_count == 4
    assert report.analysis_taxa == ["A", "B", "C", "D", "E", "F"]
    assert len(report.tree_rows) == 5
    assert len(report.coefficient_rows) == 10
    assert len(report.coefficient_summaries) == 2

    intercept = next(
        row for row in report.coefficient_summaries if row.term == "intercept"
    )
    predictor = next(
        row for row in report.coefficient_summaries if row.term == "predictor_two"
    )
    assert intercept.conclusion_stability == "stable_supported"
    assert intercept.significant_tree_count == 5
    assert predictor.dominant_direction == "positive"
    assert predictor.direction_consistency == 1.0
    assert predictor.significant_tree_count == 1
    assert predictor.significance_fraction == 0.2
    assert predictor.conclusion_stability == "mixed_support"
    assert predictor.minimum_p_value < 0.1
    assert predictor.maximum_p_value > 0.2


def test_run_posterior_tree_pgls_applies_burnin_fraction() -> None:
    report = run_posterior_tree_pgls(
        fixture("example_posterior_tree_set_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        formula="response_growth ~ predictor_two",
        lambda_value="estimate",
        burnin_fraction=0.2,
        significance_threshold=0.1,
    )

    assert report.total_tree_count == 5
    assert report.burnin_tree_count == 1
    assert report.kept_tree_count == 4
    assert report.tree_rows[0].source_tree_index == 2
    assert report.tree_rows[0].post_burnin_index == 1
    assert len(report.coefficient_rows) == 8


def test_write_posterior_tree_pgls_tables_write_expected_rows(tmp_path: Path) -> None:
    report = run_posterior_tree_pgls(
        fixture("example_posterior_tree_set_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        formula="response_growth ~ predictor_two",
        lambda_value="estimate",
        significance_threshold=0.1,
    )
    tree_out = tmp_path / "posterior-tree-pgls-trees.tsv"
    coefficient_out = tmp_path / "posterior-tree-pgls-coefficients.tsv"
    summary_out = tmp_path / "posterior-tree-pgls-summary.tsv"

    write_posterior_tree_pgls_tree_table(tree_out, report)
    write_posterior_tree_pgls_coefficient_table(coefficient_out, report)
    write_posterior_tree_pgls_summary_table(summary_out, report)

    tree_rows = tree_out.read_text(encoding="utf-8").splitlines()
    coefficient_rows = coefficient_out.read_text(encoding="utf-8").splitlines()
    summary_rows = summary_out.read_text(encoding="utf-8").splitlines()
    assert tree_rows[0].startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id"
    )
    assert coefficient_rows[0].startswith(
        "source_tree_index\tpost_burnin_index\trooted_topology_id\tterm"
    )
    assert summary_rows[0].startswith("term\ttree_fit_count\tpositive_tree_count")
    assert len(tree_rows) == 6
    assert len(coefficient_rows) == 11
    assert len(summary_rows) == 3
