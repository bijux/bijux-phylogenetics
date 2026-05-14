from __future__ import annotations

import math
from pathlib import Path
from statistics import correlation, covariance, variance

from bijux_phylogenetics.comparative.multivariate_regression import (
    MULTIVARIATE_MISSING_VALUE_POLICY,
    MULTIVARIATE_NUMERICAL_TOLERANCE,
    run_multivariate_comparative_regression,
    write_multivariate_excluded_taxa_table,
    write_multivariate_residual_association_table,
    write_multivariate_residual_correlation_table,
    write_multivariate_residual_covariance_table,
    write_multivariate_response_coefficient_table,
    write_multivariate_response_model_table,
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


def test_run_multivariate_comparative_regression_reports_covariance_and_association() -> (
    None
):
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    assert report.responses == ["response_growth", "response_range"]
    assert report.predictors == ["predictor_one", "predictor_two"]
    assert report.analysis_taxa == ["A", "B", "C", "D", "E", "F"]
    assert report.missing_value_policy == MULTIVARIATE_MISSING_VALUE_POLICY
    assert report.numerical_tolerance == MULTIVARIATE_NUMERICAL_TOLERANCE
    assert report.excluded_taxa == []
    assert len(report.response_models) == 2
    assert len(report.response_model_rows) == 2
    assert len(report.coefficient_rows) == 6
    assert len(report.covariance_rows) == 4
    assert len(report.correlation_rows) == 4
    assert len(report.association_rows) == 1
    association = report.association_rows[0]
    assert association.left_response == "response_growth"
    assert association.right_response == "response_range"
    assert association.pair_count == 6
    assert math.isclose(association.correlation, 0.0, abs_tol=1e-12)
    assert association.p_value == 1.0
    assert report.covariance_diagnostics.response_count == 2
    assert report.covariance_diagnostics.matrix_rank == 1
    assert report.covariance_diagnostics.is_singular is True
    assert report.covariance_diagnostics.is_near_singular is True
    assert math.isinf(report.covariance_diagnostics.condition_number)
    diagonal = next(
        row
        for row in report.covariance_rows
        if row.left_response == "response_growth"
        and row.right_response == "response_growth"
    )
    assert diagonal.is_diagonal is True
    assert math.isclose(diagonal.correlation, 1.0)
    assert math.isfinite(report.response_model_rows[0].log_likelihood)
    assert report.coefficient_rows[0].standard_error >= 0.0


def test_run_multivariate_comparative_regression_matches_independent_python_reference() -> (
    None
):
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multiple.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    left_model = report.response_models[0]
    right_model = report.response_models[1]
    diagonal = next(
        row
        for row in report.covariance_rows
        if row.left_response == left_model.response
        and row.right_response == left_model.response
    )
    pair = report.association_rows[0]
    assert math.isclose(
        diagonal.covariance,
        variance(left_model.residuals),
        abs_tol=report.numerical_tolerance,
    )
    assert math.isclose(
        pair.covariance,
        covariance(left_model.residuals, right_model.residuals),
        abs_tol=report.numerical_tolerance,
    )
    assert math.isclose(
        pair.correlation,
        correlation(left_model.residuals, right_model.residuals),
        abs_tol=report.numerical_tolerance,
    )


def test_run_multivariate_comparative_regression_reports_full_rank_covariance_diagnostics() -> (
    None
):
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multivariate_full_rank.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    assert report.covariance_diagnostics.response_count == 2
    assert report.covariance_diagnostics.matrix_rank == 2
    assert report.covariance_diagnostics.is_singular is False
    assert report.covariance_diagnostics.is_near_singular is False
    assert math.isfinite(report.covariance_diagnostics.condition_number)


def test_run_multivariate_comparative_regression_excludes_incomplete_taxa() -> None:
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multivariate_missing.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    assert report.analysis_taxa == ["A", "C", "D", "E", "F"]
    assert report.missing_value_policy == MULTIVARIATE_MISSING_VALUE_POLICY
    assert len(report.excluded_taxa) == 1
    excluded = report.excluded_taxa[0]
    assert excluded.taxon == "B"
    assert excluded.reason == "missing_required_values"
    assert excluded.missing_columns == ["response_range"]
    assert excluded.blocking_responses == ["response_range"]


def test_run_multivariate_comparative_regression_supports_categorical_and_interaction_terms() -> (
    None
):
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multivariate_interaction.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "habitat", "predictor_one:habitat"],
        lambda_value=0.0,
    )
    assert report.analysis_taxa == ["A", "B", "C", "D", "E", "F"]
    assert all(
        model.formula.interaction_terms == ["predictor_one:habitat"]
        for model in report.response_models
    )
    encoded_terms = {row.term for row in report.coefficient_rows}
    assert "habitat[tundra]" in encoded_terms
    assert "predictor_one:habitat[tundra]" in encoded_terms
    assert any(
        "weak residual degrees of freedom" in warning for warning in report.warnings
    )


def test_run_multivariate_comparative_regression_warns_on_singular_covariance() -> None:
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multivariate_singular.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one"],
        lambda_value=0.0,
    )
    assert any("singular or near-singular" in warning for warning in report.warnings)
    diagonal_rows = [row for row in report.covariance_rows if row.is_diagonal]
    assert all(
        math.isclose(row.covariance, 0.0, abs_tol=report.numerical_tolerance)
        for row in diagonal_rows
    )
    assert report.covariance_diagnostics.matrix_rank == 0
    assert report.covariance_diagnostics.is_singular is True
    assert report.covariance_diagnostics.is_near_singular is True
    assert math.isinf(report.covariance_diagnostics.condition_number)


def test_run_multivariate_comparative_regression_detects_collinear_residual_covariance() -> (
    None
):
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multivariate_collinear.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one"],
        lambda_value=0.0,
    )
    diagonal_rows = [row for row in report.covariance_rows if row.is_diagonal]
    assert all(row.covariance > report.numerical_tolerance for row in diagonal_rows)
    assert any("singular or near-singular" in warning for warning in report.warnings)
    assert report.covariance_diagnostics.matrix_rank == 1
    assert report.covariance_diagnostics.is_singular is True
    assert report.covariance_diagnostics.is_near_singular is True
    assert math.isinf(report.covariance_diagnostics.condition_number)


def test_write_multivariate_regression_tables_write_review_ledgers(
    tmp_path: Path,
) -> None:
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multivariate_missing.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    model_path = tmp_path / "multivariate-response-models.tsv"
    coefficient_path = tmp_path / "multivariate-response-coefficients.tsv"
    covariance_path = tmp_path / "multivariate-residual-covariance.tsv"
    correlation_path = tmp_path / "multivariate-residual-correlation.tsv"
    association_path = tmp_path / "multivariate-residual-associations.tsv"
    excluded_path = tmp_path / "multivariate-excluded-taxa.tsv"
    write_multivariate_response_model_table(model_path, report)
    write_multivariate_response_coefficient_table(coefficient_path, report)
    write_multivariate_residual_covariance_table(covariance_path, report)
    write_multivariate_residual_correlation_table(correlation_path, report)
    write_multivariate_residual_association_table(association_path, report)
    write_multivariate_excluded_taxa_table(excluded_path, report)
    model_rows = model_path.read_text(encoding="utf-8").splitlines()
    coefficient_rows = coefficient_path.read_text(encoding="utf-8").splitlines()
    covariance_rows = covariance_path.read_text(encoding="utf-8").splitlines()
    correlation_rows = correlation_path.read_text(encoding="utf-8").splitlines()
    association_rows = association_path.read_text(encoding="utf-8").splitlines()
    excluded_rows = excluded_path.read_text(encoding="utf-8").splitlines()
    assert model_rows[0].startswith("response\tformula\tpredictor_term_count")
    assert coefficient_rows[0].startswith("response\tformula\tterm\testimate")
    assert covariance_rows[0].startswith(
        "left_response\tright_response\tpair_count\tis_diagonal"
    )
    assert correlation_rows[0].startswith(
        "left_response\tright_response\tpair_count\tis_diagonal\tcorrelation"
    )
    assert association_rows[0].startswith(
        "left_response\tright_response\tpair_count\tcovariance\tcorrelation"
    )
    assert excluded_rows[0] == (
        "taxon\treason\tmissing_columns\tblocking_responses\tdetails"
    )
    assert len(model_rows) == 3
    assert len(coefficient_rows) == 7
    assert len(covariance_rows) == 5
    assert len(correlation_rows) == 5
    assert len(association_rows) == 2
    assert len(excluded_rows) == 2
