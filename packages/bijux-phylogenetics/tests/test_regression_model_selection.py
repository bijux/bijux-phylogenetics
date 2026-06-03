from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.comparative.regression import (
    compare_comparative_regression_models,
    write_comparative_regression_excluded_taxa_table,
    write_comparative_regression_model_ranking_table,
    write_comparative_regression_pairwise_table,
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


def test_compare_comparative_regression_models_ranks_continuous_candidates() -> None:
    report = compare_comparative_regression_models(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        formulas=[
            "response_growth ~ predictor_one",
            "response_growth ~ predictor_two",
            "response_growth ~ predictor_one + predictor_two",
        ],
        lambda_value=0.0,
    )
    assert report.model_family == "pgls"
    assert report.selected_criterion == "AICc"
    assert report.analysis_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.best_formula == "response_growth ~ predictor_one"
    assert [row.formula for row in sorted(report.rows, key=lambda row: row.rank)] == [
        "response_growth ~ predictor_one",
        "response_growth ~ predictor_one + predictor_two",
        "response_growth ~ predictor_two",
    ]
    pairwise = {
        (row.left_formula, row.right_formula): row.comparison_kind
        for row in report.pairwise_rows
    }
    assert (
        pairwise[
            (
                "response_growth ~ predictor_one",
                "response_growth ~ predictor_one + predictor_two",
            )
        ]
        == "left_nested_in_right"
    )
    assert (
        pairwise[
            (
                "response_growth ~ predictor_one",
                "response_growth ~ predictor_two",
            )
        ]
        == "non_nested"
    )


def test_compare_comparative_regression_models_matches_logistic_reference_ranking() -> (
    None
):
    reference = json.loads(
        fixture("comparative_regression_model_selection_reference.json").read_text(
            encoding="utf-8"
        )
    )
    report = compare_comparative_regression_models(
        fixture("example_tree_eight_taxa.nwk"),
        fixture("example_traits_phylogenetic_logistic_model_selection.tsv"),
        formulas=[
            "presence ~ body_size",
            "presence ~ habitat",
            "presence ~ body_size + habitat",
        ],
        lambda_value=1.0,
    )
    ranked_formulas = [
        row.formula for row in sorted(report.rows, key=lambda row: row.rank)
    ]
    reference_ranked = [
        row["formula"]
        for row in sorted(reference["rows"], key=lambda row: (row["aic"], row["bic"]))
    ]
    assert report.model_family == "logistic"
    assert report.best_formula == reference["best_formula"]
    assert ranked_formulas == reference_ranked
    assert any(row.separation_detected for row in report.rows)
    assert "selected model" in " ".join(report.warnings)


def test_write_comparative_regression_model_tables_write_rows(tmp_path: Path) -> None:
    report = compare_comparative_regression_models(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        formulas=[
            "response_growth ~ predictor_one",
            "response_growth ~ predictor_two",
        ],
        lambda_value=0.0,
    )
    ranking_out = tmp_path / "comparative-model-ranking.tsv"
    pairwise_out = tmp_path / "comparative-model-pairwise.tsv"
    excluded_out = tmp_path / "comparative-model-excluded.tsv"
    write_comparative_regression_model_ranking_table(ranking_out, report)
    write_comparative_regression_pairwise_table(pairwise_out, report)
    write_comparative_regression_excluded_taxa_table(excluded_out, report)
    ranking_rows = ranking_out.read_text(encoding="utf-8").splitlines()
    pairwise_rows = pairwise_out.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_out.read_text(encoding="utf-8").splitlines()
    assert ranking_rows[0].startswith("formula\tmodel_family\tparameter_count")
    assert len(ranking_rows) == 3
    assert pairwise_rows[0].startswith("left_formula\tright_formula\tcomparison_kind")
    assert len(pairwise_rows) == 2
    assert excluded_rows == ["taxon\treason\tmissing_columns"]


def test_compare_comparative_regression_models_requires_shared_response() -> None:
    try:
        compare_comparative_regression_models(
            fixture("example_tree_six_taxa.nwk"),
            fixture("example_traits_comparative_multiple.tsv"),
            formulas=[
                "response_growth ~ predictor_one",
                "response_range ~ predictor_one",
            ],
            lambda_value=0.0,
        )
    except Exception as error:  # pragma: no cover - tightened by assertion below
        message = str(error)
    else:  # pragma: no cover
        raise AssertionError("expected shared-response validation to fail")
    assert "same response column" in message
