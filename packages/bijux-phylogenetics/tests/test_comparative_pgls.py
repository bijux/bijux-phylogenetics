from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import (
    inspect_pgls_inputs,
    run_pgls,
    run_pgls_multiple_testing,
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


def test_pgls_input_inspection_encodes_categorical_predictors() -> None:
    report = inspect_pgls_inputs(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one", "habitat"],
    )
    assert report.ready is True
    assert report.categorical_predictors == ["habitat"]
    habitat = next(row for row in report.predictors if row.name == "habitat")
    assert habitat.reference_level == "forest"
    assert habitat.encoded_columns == ["habitat[tundra]"]
    assert habitat.observed_levels == ["forest", "tundra"]
    assert habitat.level_counts == {"forest": 2, "tundra": 2}
    assert report.encoded_columns == ["intercept", "predictor_one", "habitat[tundra]"]
    assert report.formula_audit.parameter_count == 3
    assert report.formula_audit.minimum_required_taxa == 4


def test_pgls_input_inspection_reports_transformed_terms_and_exclusions() -> None:
    report = inspect_pgls_inputs(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        formula="response ~ log(predictor_one) + habitat",
    )
    transformed = next(
        row for row in report.predictors if row.name == "log(predictor_one)"
    )
    assert transformed.transformation == "log"
    assert transformed.source_column == "predictor_one"
    assert report.formula_audit.transformed_terms == ["log(predictor_one)"]


def test_pgls_input_inspection_rejects_missing_branch_lengths() -> None:
    report = inspect_pgls_inputs(
        fixture("example_tree_no_lengths.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
    )
    assert report.ready is False
    assert "PGLS requires complete tree branch lengths" in report.blockers


def test_run_pgls_supports_one_predictor() -> None:
    report = run_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one"],
        lambda_value=0.0,
    )
    assert report.predictors == ["predictor_one"]
    assert report.taxon_count == 4
    assert len(report.coefficients) == 2
    assert math.isclose(report.coefficients[0].estimate, 1.0, abs_tol=1e-6)
    assert math.isclose(report.coefficients[1].estimate, 0.7, abs_tol=1e-6)


def test_run_pgls_supports_multiple_predictors() -> None:
    report = run_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    coefficients = {
        coefficient.name: coefficient.estimate for coefficient in report.coefficients
    }
    assert math.isclose(coefficients["intercept"], 1.0, abs_tol=1e-6)
    assert math.isclose(coefficients["predictor_one"], 0.5, abs_tol=1e-6)
    assert math.isclose(coefficients["predictor_two"], 1.0, abs_tol=1e-6)


def test_run_pgls_supports_categorical_predictors_with_dummy_encoding() -> None:
    report = run_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        response="response",
        predictors=["predictor_one", "habitat"],
        lambda_value=0.0,
    )
    coefficients = {
        coefficient.name: coefficient.estimate for coefficient in report.coefficients
    }
    assert report.encoded_columns == ["intercept", "predictor_one", "habitat[tundra]"]
    assert math.isclose(coefficients["predictor_one"], 1.5, abs_tol=1e-6)
    assert math.isclose(coefficients["habitat[tundra]"], -2.0, abs_tol=1e-6)
    assert report.diagnostics.outlier_taxa == []


def test_pgls_formula_expands_interactions() -> None:
    report = inspect_pgls_inputs(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_interaction.tsv"),
        formula="response ~ predictor_one * habitat",
    )
    assert report.response == "response"
    assert report.formula.predictors == ["predictor_one", "habitat"]
    assert report.formula.interaction_terms == ["predictor_one:habitat"]
    assert report.encoded_columns == [
        "intercept",
        "predictor_one",
        "habitat[tundra]",
        "predictor_one:habitat[tundra]",
    ]


def test_run_pgls_supports_formula_interactions() -> None:
    report = run_pgls(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_interaction.tsv"),
        formula="response ~ predictor_one * habitat",
        lambda_value=0.0,
    )
    coefficients = {
        coefficient.name: coefficient.estimate for coefficient in report.coefficients
    }
    assert report.formula.formula == "response ~ predictor_one * habitat"
    assert math.isclose(coefficients["intercept"], 1.0, abs_tol=1e-6)
    assert math.isclose(coefficients["predictor_one"], 1.0, abs_tol=1e-6)
    assert math.isclose(coefficients["habitat[tundra]"], 2.0, abs_tol=1e-6)
    assert math.isclose(
        coefficients["predictor_one:habitat[tundra]"], 0.5, abs_tol=1e-6
    )


def test_run_pgls_supports_transformed_numeric_predictors() -> None:
    report = run_pgls(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        formula="response ~ log(predictor_one) + habitat",
        lambda_value=0.0,
    )
    coefficients = {
        coefficient.name: coefficient.estimate for coefficient in report.coefficients
    }
    assert "log(predictor_one)" in coefficients
    assert "habitat[tundra]" in coefficients


def test_pgls_overfit_guard_blocks_saturated_interaction_model() -> None:
    report = inspect_pgls_inputs(
        fixture("example_tree.nwk"),
        fixture("example_traits_comparative.tsv"),
        formula="response ~ predictor_one * habitat",
    )
    assert report.ready is False
    assert (
        "PGLS overfit guard requires at least one residual degree of freedom after predictor encoding"
        in report.blockers
    )


def test_run_pgls_multiple_testing_adjusts_p_values() -> None:
    report = run_pgls_multiple_testing(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    assert report.adjustment_method == "benjamini-hochberg"
    assert report.family_size == 4
    assert len(report.rows) == 4
    assert report.raw_significant_count >= report.adjusted_significant_count
    assert all(row.adjusted_p_value >= row.p_value for row in report.rows)
    assert all(0.0 <= row.adjusted_p_value <= 1.0 for row in report.rows)
