from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy

from bijux_phylogenetics.comparative._math import student_t_two_sided_p_value
from bijux_phylogenetics.comparative.common import (
    build_brownian_covariance_matrix,
    lambda_transform_covariance,
)
from bijux_phylogenetics.comparative.regression import (
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
from bijux_phylogenetics.io.trees import load_tree

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")
REFERENCE_PARITY_TOLERANCE = 2e-9


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def _load_numeric_rows(path: Path) -> dict[str, dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {
            row["taxon"]: {
                key: float(value)
                for key, value in row.items()
                if key != "taxon" and value is not None and value != ""
            }
            for row in reader
        }


def _multivariate_reference_gls(
    *,
    tree_path: Path,
    traits_path: Path,
    responses: list[str],
    predictors: list[str],
    taxa: list[str],
    lambda_value: float | dict[str, float],
) -> dict[str, object]:
    rows_by_taxon = _load_numeric_rows(traits_path)
    tree = load_tree(tree_path)
    base_covariance = build_brownian_covariance_matrix(tree, taxa)
    design_matrix = numpy.array(
        [
            [1.0, *(rows_by_taxon[taxon][predictor] for predictor in predictors)]
            for taxon in taxa
        ],
        dtype=float,
    )
    coefficient_names = ["intercept", *predictors]
    coefficients_by_response: dict[str, dict[str, tuple[float, float, float]]] = {}
    residual_matrix_rows: list[list[float]] = []
    log_likelihoods: dict[str, float] = {}
    residual_variances: dict[str, float] = {}
    for response in responses:
        resolved_lambda = (
            lambda_value[response] if isinstance(lambda_value, dict) else lambda_value
        )
        covariance = numpy.array(
            lambda_transform_covariance(base_covariance, resolved_lambda),
            dtype=float,
        )
        inverse_covariance = numpy.linalg.inv(covariance)
        sign, log_determinant = numpy.linalg.slogdet(covariance)
        assert sign > 0.0
        xt_vinv_x = design_matrix.T @ inverse_covariance @ design_matrix
        covariance_of_betas = numpy.linalg.inv(xt_vinv_x)
        observed = numpy.array(
            [rows_by_taxon[taxon][response] for taxon in taxa],
            dtype=float,
        )
        coefficients = covariance_of_betas @ (
            design_matrix.T @ inverse_covariance @ observed
        )
        fitted = design_matrix @ coefficients
        residuals = observed - fitted
        residual_matrix_rows.append(residuals.tolist())
        degrees_of_freedom = len(taxa) - len(coefficient_names)
        residual_scale = float(residuals.T @ inverse_covariance @ residuals)
        residual_variance = residual_scale / degrees_of_freedom
        residual_variances[response] = residual_variance
        log_likelihoods[response] = -0.5 * (
            len(taxa) * math.log(2.0 * math.pi * max(residual_scale / len(taxa), 1e-12))
            + float(log_determinant)
            + len(taxa)
        )
        coefficients_by_response[response] = {}
        for index, term in enumerate(coefficient_names):
            standard_error = math.sqrt(
                max(float(covariance_of_betas[index, index]) * residual_variance, 0.0)
            )
            test_statistic = (
                float(coefficients[index]) / standard_error if standard_error else 0.0
            )
            coefficients_by_response[response][term] = (
                float(coefficients[index]),
                standard_error,
                student_t_two_sided_p_value(test_statistic, degrees_of_freedom),
            )
    residual_matrix = numpy.array(residual_matrix_rows, dtype=float)
    return {
        "coefficients_by_response": coefficients_by_response,
        "log_likelihoods": log_likelihoods,
        "residual_variances": residual_variances,
        "residual_covariance": numpy.cov(residual_matrix, ddof=1),
        "residual_correlation": numpy.corrcoef(residual_matrix),
    }


def _assert_multivariate_reference_matches_report(
    report, reference: dict[str, object]
) -> None:
    for row in report.response_model_rows:
        assert math.isclose(
            row.log_likelihood,
            reference["log_likelihoods"][row.response],
            abs_tol=REFERENCE_PARITY_TOLERANCE,
        )
        assert math.isclose(
            row.residual_variance,
            reference["residual_variances"][row.response],
            abs_tol=REFERENCE_PARITY_TOLERANCE,
        )
    coefficient_reference = reference["coefficients_by_response"]
    for row in report.coefficient_rows:
        estimate, standard_error, p_value = coefficient_reference[row.response][
            row.term
        ]
        assert math.isclose(
            row.estimate,
            estimate,
            abs_tol=REFERENCE_PARITY_TOLERANCE,
        )
        assert math.isclose(
            row.standard_error,
            standard_error,
            abs_tol=REFERENCE_PARITY_TOLERANCE,
        )
        assert math.isclose(row.p_value, p_value, abs_tol=REFERENCE_PARITY_TOLERANCE)
    response_index = {
        response: index for index, response in enumerate(report.responses)
    }
    covariance_reference = reference["residual_covariance"]
    correlation_reference = reference["residual_correlation"]
    for row in report.covariance_rows:
        left_index = response_index[row.left_response]
        right_index = response_index[row.right_response]
        assert math.isclose(
            row.covariance,
            float(covariance_reference[left_index, right_index]),
            abs_tol=REFERENCE_PARITY_TOLERANCE,
        )
        assert math.isclose(
            row.correlation,
            float(correlation_reference[left_index, right_index]),
            abs_tol=REFERENCE_PARITY_TOLERANCE,
        )
    for row in report.correlation_rows:
        left_index = response_index[row.left_response]
        right_index = response_index[row.right_response]
        assert math.isclose(
            row.correlation,
            float(correlation_reference[left_index, right_index]),
            abs_tol=REFERENCE_PARITY_TOLERANCE,
        )


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
        fixture("example_traits_comparative_multivariate_full_rank.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.0,
    )
    reference = _multivariate_reference_gls(
        tree_path=fixture("example_tree_six_taxa.nwk"),
        traits_path=fixture("example_traits_comparative_multivariate_full_rank.tsv"),
        responses=report.responses,
        predictors=report.predictors,
        taxa=report.analysis_taxa,
        lambda_value=0.0,
    )
    _assert_multivariate_reference_matches_report(report, reference)


def test_run_multivariate_comparative_regression_matches_independent_python_reference_with_nonzero_lambda() -> (
    None
):
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multivariate_full_rank.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value=0.5,
    )
    reference = _multivariate_reference_gls(
        tree_path=fixture("example_tree_six_taxa.nwk"),
        traits_path=fixture("example_traits_comparative_multivariate_full_rank.tsv"),
        responses=report.responses,
        predictors=report.predictors,
        taxa=report.analysis_taxa,
        lambda_value=0.5,
    )
    _assert_multivariate_reference_matches_report(report, reference)


def test_run_multivariate_comparative_regression_warns_when_estimated_lambdas_materially_diverge() -> (
    None
):
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multivariate_heterogeneous_lambda.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value="estimate",
    )
    response_lambda = {
        row.response: row.lambda_value for row in report.response_model_rows
    }
    assert response_lambda == {"response_growth": 1.0, "response_range": 0.0}
    assert any(
        "materially different Pagel lambda values (0 to 1)" in warning
        for warning in report.warnings
    )


def test_run_multivariate_comparative_regression_matches_reference_with_heterogeneous_estimated_lambdas() -> (
    None
):
    report = run_multivariate_comparative_regression(
        fixture("example_tree_six_taxa.nwk"),
        fixture("example_traits_comparative_multivariate_heterogeneous_lambda.tsv"),
        responses=["response_growth", "response_range"],
        predictors=["predictor_one", "predictor_two"],
        lambda_value="estimate",
    )
    reference = _multivariate_reference_gls(
        tree_path=fixture("example_tree_six_taxa.nwk"),
        traits_path=fixture(
            "example_traits_comparative_multivariate_heterogeneous_lambda.tsv"
        ),
        responses=report.responses,
        predictors=report.predictors,
        taxa=report.analysis_taxa,
        lambda_value={
            model.response: model.lambda_value for model in report.response_models
        },
    )
    _assert_multivariate_reference_matches_report(report, reference)


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
