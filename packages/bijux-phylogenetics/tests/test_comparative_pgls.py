from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

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
    intercept = next(row for row in report.coefficients if row.name == "intercept")
    slope = next(row for row in report.coefficients if row.name == "predictor_one")
    assert math.isclose(intercept.estimate, 1.0, abs_tol=1e-6)
    assert math.isclose(slope.estimate, 0.7, abs_tol=1e-6)
    assert math.isclose(
        report.aic,
        -2.0 * report.log_likelihood + (2.0 * (len(report.coefficients) + 1)),
        rel_tol=1e-12,
    )
    assert all(row.inference_distribution == "student-t" for row in report.coefficients)
    assert all(row.degrees_of_freedom == 2 for row in report.coefficients)
    assert math.isclose(slope.test_statistic, 2.4748737341529177, rel_tol=1e-12)
    assert math.isclose(slope.p_value, 0.1317568578755406, rel_tol=1e-12)
    assert math.isclose(
        slope.lower_95_confidence_interval, -0.5169739689186594, rel_tol=1e-12
    )
    assert math.isclose(
        slope.upper_95_confidence_interval, 1.9169739689186598, rel_tol=1e-12
    )
    assert slope.lower_95_confidence_interval < slope.estimate
    assert slope.upper_95_confidence_interval > slope.estimate


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


def test_run_pgls_matches_primate_reference_when_lambda_is_fixed() -> None:
    repository_root = Path(__file__).resolve().parents[3]
    parity_fixture = (
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-pgls-and-signal"
        / "evidence-003"
        / "results"
        / "pagel-lambda-regression-parity.json"
    )
    reference = json.loads(parity_fixture.read_text(encoding="utf-8"))
    lambda_value = reference["r_fixed_reference_lambda"]["lambda_value"]
    report = run_pgls(
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_trimmed_primatetree.nwk",
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_primate.csv",
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value=lambda_value,
    )
    coefficients = {row.name: row for row in report.coefficients}
    fixed_reference = reference["r_fixed_reference_lambda"]
    assert math.isclose(report.lambda_value, lambda_value, abs_tol=1e-12)
    assert math.isclose(
        report.log_likelihood,
        fixed_reference["log_likelihood"],
        rel_tol=5e-4,
        abs_tol=5e-4,
    )
    assert math.isclose(
        report.aic,
        fixed_reference["aic"],
        rel_tol=5e-4,
        abs_tol=5e-4,
    )
    assert math.isclose(
        coefficients["intercept"].estimate,
        fixed_reference["coefficients"]["intercept"],
        rel_tol=5e-4,
        abs_tol=5e-4,
    )
    assert math.isclose(
        coefficients["social_group_size"].estimate,
        fixed_reference["coefficients"]["social_group_size"],
        rel_tol=5e-4,
        abs_tol=5e-4,
    )
    assert math.isclose(
        coefficients["intercept"].standard_error,
        fixed_reference["standard_errors"]["intercept"],
        rel_tol=5e-4,
        abs_tol=5e-4,
    )
    assert math.isclose(
        coefficients["social_group_size"].standard_error,
        fixed_reference["standard_errors"]["social_group_size"],
        rel_tol=5e-4,
        abs_tol=5e-4,
    )
    assert math.isclose(
        coefficients["intercept"].p_value,
        fixed_reference["p_values"]["intercept"],
        rel_tol=5e-4,
        abs_tol=5e-7,
    )
    assert math.isclose(
        coefficients["social_group_size"].p_value,
        fixed_reference["p_values"]["social_group_size"],
        rel_tol=5e-4,
        abs_tol=5e-6,
    )


@pytest.mark.slow
def test_run_pgls_reports_aic_for_estimated_lambda() -> None:
    repository_root = Path(__file__).resolve().parents[3]
    report = run_pgls(
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_trimmed_primatetree.nwk",
        repository_root
        / "evidence-book"
        / "studies"
        / "primate-longevity-signal"
        / "datasets"
        / "reference_primate.csv",
        response="longevity",
        predictors=["social_group_size"],
        taxon_column="species",
        lambda_value="estimate",
    )
    assert math.isclose(
        report.aic,
        -2.0 * report.log_likelihood + (2.0 * (len(report.coefficients) + 2)),
        rel_tol=1e-12,
    )
