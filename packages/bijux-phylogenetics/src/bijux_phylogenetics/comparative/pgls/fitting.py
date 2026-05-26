from __future__ import annotations

from collections.abc import Callable
import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    matrix_multiply,
    student_t_quantile,
    student_t_two_sided_p_value,
    transpose,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    lambda_transform_covariance,
    load_comparative_dataset,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .design import inspect_pgls_inputs
from .models import (
    PGLSCoefficient,
    PGLSDiagnosticsReport,
    PGLSFittedObservation,
    PGLSLambdaFitReport,
    PGLSLambdaProfileRow,
    PGLSLeverageRow,
    PGLSResidualOutlier,
    PGLSResult,
)


def run_pgls(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> PGLSResult:
    """Fit a phylogenetic generalized least-squares model."""
    input_report = inspect_pgls_inputs(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
    )
    if not input_report.ready:
        raise ComparativeMethodError("; ".join(input_report.blockers))

    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=input_report.formula_audit.response_column,
        taxon_column=taxon_column,
        minimum_taxa=len(input_report.encoded_columns) + 1,
        require_rooted=True,
        require_binary=False,
    )
    taxa = list(input_report.analysis_taxa)
    response_values = [row.response_value for row in input_report.model_matrix.rows]
    encoded_columns = list(input_report.model_matrix.encoded_columns)
    design_matrix = [
        [row.encoded_values[column] for column in encoded_columns]
        for row in input_report.model_matrix.rows
    ]
    dataset = ComparativeDataset(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        tree=dataset.tree,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        taxa=taxa,
        trait_values=response_values,
        covariance_matrix=_subset_covariance(
            dataset.covariance_matrix, dataset.taxa, taxa
        ),
        readiness=dataset.readiness,
    )
    lambda_fit = _resolve_lambda_fit(
        dataset,
        design_matrix,
        response_values,
        lambda_value,
    )
    resolved_lambda = lambda_fit.lambda_value
    covariance = lambda_transform_covariance(dataset.covariance_matrix, resolved_lambda)
    inverse_covariance = invert_matrix(covariance)
    coefficients, covariance_of_betas, fitted_values = _fit_gls(
        design_matrix, response_values, inverse_covariance
    )
    residuals = [
        observed - fitted
        for observed, fitted in zip(response_values, fitted_values, strict=True)
    ]
    degrees_of_freedom = len(response_values) - len(coefficients)
    residual_variance = (
        _quadratic_form(residuals, inverse_covariance) / degrees_of_freedom
    )
    coefficient_reports: list[PGLSCoefficient] = []
    critical_value = student_t_quantile(0.975, degrees_of_freedom)
    for index, name in enumerate(encoded_columns):
        standard_error = math.sqrt(
            max(covariance_of_betas[index][index] * residual_variance, 0.0)
        )
        test_statistic = coefficients[index] / standard_error if standard_error else 0.0
        p_value = student_t_two_sided_p_value(test_statistic, degrees_of_freedom)
        interval_radius = critical_value * standard_error
        coefficient_reports.append(
            PGLSCoefficient(
                name=name,
                estimate=coefficients[index],
                standard_error=standard_error,
                test_statistic=test_statistic,
                p_value=p_value,
                lower_95_confidence_interval=coefficients[index] - interval_radius,
                upper_95_confidence_interval=coefficients[index] + interval_radius,
                degrees_of_freedom=degrees_of_freedom,
                inference_distribution="student-t",
            )
        )
    mean_response = sum(response_values) / len(response_values)
    total_sum_of_squares = sum(
        (value - mean_response) ** 2 for value in response_values
    )
    residual_sum_of_squares = sum(value * value for value in residuals)
    r_squared = (
        1.0 - (residual_sum_of_squares / total_sum_of_squares)
        if total_sum_of_squares
        else 1.0
    )
    log_likelihood = _gls_log_likelihood(
        response_values,
        residuals,
        inverse_covariance,
        covariance,
    )
    parameter_count = len(coefficients) + 1
    if lambda_fit.mode == "estimated":
        parameter_count += 1
    aic = -2.0 * log_likelihood + (2.0 * parameter_count)
    return PGLSResult(
        tree_path=tree_path,
        traits_path=traits_path,
        response=input_report.response,
        formula=input_report.formula,
        predictors=list(input_report.formula.predictors),
        interaction_terms=list(input_report.formula.interaction_terms),
        encoded_columns=encoded_columns,
        taxon_count=len(taxa),
        lambda_value=resolved_lambda,
        lambda_fit=lambda_fit,
        log_likelihood=log_likelihood,
        aic=aic,
        residual_variance=residual_variance,
        r_squared=r_squared,
        coefficients=coefficient_reports,
        fitted_values=fitted_values,
        residuals=residuals,
        taxa=taxa,
        diagnostics=_build_pgls_diagnostics(
            taxa,
            response_values,
            fitted_values,
            residuals,
            residual_variance,
            design_matrix,
            inverse_covariance,
        ),
    )


def _fit_gls(
    design_matrix: list[list[float]],
    response_values: list[float],
    inverse_covariance: list[list[float]],
) -> tuple[list[float], list[list[float]], list[float]]:
    x_transposed = transpose(design_matrix)
    xt_vinv = matrix_multiply(x_transposed, inverse_covariance)
    xt_vinv_x = matrix_multiply(xt_vinv, design_matrix)
    xt_vinv_x_inverse = invert_matrix(xt_vinv_x)
    xt_vinv_y = [
        sum(
            row[index] * response_values[index] for index in range(len(response_values))
        )
        for row in xt_vinv
    ]
    coefficients = [
        sum(row[index] * xt_vinv_y[index] for index in range(len(xt_vinv_y)))
        for row in xt_vinv_x_inverse
    ]
    fitted_values = [
        sum(beta * value for beta, value in zip(coefficients, row, strict=True))
        for row in design_matrix
    ]
    return coefficients, xt_vinv_x_inverse, fitted_values


def _build_pgls_diagnostics(
    taxa: list[str],
    observed_values: list[float],
    fitted_values: list[float],
    residuals: list[float],
    residual_variance: float,
    design_matrix: list[list[float]],
    inverse_covariance: list[list[float]],
) -> PGLSDiagnosticsReport:
    x_transposed = transpose(design_matrix)
    xt_vinv = matrix_multiply(x_transposed, inverse_covariance)
    xt_vinv_x_inverse = invert_matrix(matrix_multiply(xt_vinv, design_matrix))
    hat = matrix_multiply(design_matrix, matrix_multiply(xt_vinv_x_inverse, xt_vinv))
    leverage_rows: list[PGLSLeverageRow] = []
    outlier_taxa: list[PGLSResidualOutlier] = []
    fitted_observed_rows: list[PGLSFittedObservation] = []
    residual_mean = sum(residuals) / len(residuals)
    for index, taxon in enumerate(taxa):
        leverage = min(max(hat[index][index], 0.0), 0.999999)
        denominator = math.sqrt(max(residual_variance * (1.0 - leverage), 1e-12))
        standardized = residuals[index] / denominator
        leverage_rows.append(
            PGLSLeverageRow(
                taxon=taxon,
                leverage=leverage,
                standardized_residual=standardized,
            )
        )
        fitted_observed_rows.append(
            PGLSFittedObservation(
                taxon=taxon,
                observed=observed_values[index],
                fitted=fitted_values[index],
                residual=residuals[index],
            )
        )
        if abs(standardized) >= 2.0:
            outlier_taxa.append(
                PGLSResidualOutlier(
                    taxon=taxon,
                    residual=residuals[index],
                    standardized_residual=standardized,
                )
            )
    return PGLSDiagnosticsReport(
        residual_mean=residual_mean,
        leverage_rows=leverage_rows,
        outlier_taxa=outlier_taxa,
        fitted_observed_rows=fitted_observed_rows,
    )


def _resolve_lambda_fit(
    dataset: ComparativeDataset,
    design_matrix: list[list[float]],
    response_values: list[float],
    lambda_value: float | str,
) -> PGLSLambdaFitReport:
    likelihood_cache: dict[float, float] = {}

    def _cached_log_likelihood(candidate: float) -> float:
        if candidate not in likelihood_cache:
            likelihood_cache[candidate] = _lambda_log_likelihood(
                dataset,
                design_matrix,
                response_values,
                candidate,
            )
        return likelihood_cache[candidate]

    null_log_likelihood = _cached_log_likelihood(0.0)
    brownian_log_likelihood = _cached_log_likelihood(1.0)
    if isinstance(lambda_value, (float, int)):
        if not 0.0 <= lambda_value <= 1.0:
            raise ComparativeMethodError(
                "PGLS lambda must be between 0 and 1 inclusive"
            )
        resolved_lambda = float(lambda_value)
        log_likelihood = _cached_log_likelihood(resolved_lambda)
        return PGLSLambdaFitReport(
            mode="fixed",
            lambda_value=resolved_lambda,
            log_likelihood=log_likelihood,
            null_log_likelihood=null_log_likelihood,
            brownian_log_likelihood=brownian_log_likelihood,
            lower_95_confidence_interval=None,
            upper_95_confidence_interval=None,
            profile_rows=[
                PGLSLambdaProfileRow(
                    lambda_value=resolved_lambda,
                    log_likelihood=log_likelihood,
                    delta_log_likelihood=0.0,
                    within_95_confidence_interval=True,
                )
            ],
        )
    if lambda_value != "estimate":
        raise ComparativeMethodError(
            "PGLS lambda must be 'estimate' or a numeric value"
        )
    return _estimate_lambda_for_pgls(
        log_likelihood_at_lambda=_cached_log_likelihood,
        null_log_likelihood=null_log_likelihood,
        brownian_log_likelihood=brownian_log_likelihood,
    )


def _subset_covariance(
    covariance_matrix: list[list[float]],
    source_taxa: list[str],
    target_taxa: list[str],
) -> list[list[float]]:
    positions = {taxon: index for index, taxon in enumerate(source_taxa)}
    return [
        [
            covariance_matrix[positions[left_taxon]][positions[right_taxon]]
            for right_taxon in target_taxa
        ]
        for left_taxon in target_taxa
    ]


def _quadratic_form(vector: list[float], matrix: list[list[float]]) -> float:
    total = 0.0
    for row_index, row in enumerate(matrix):
        total += vector[row_index] * sum(
            value * vector[column_index] for column_index, value in enumerate(row)
        )
    return total


def _gls_log_likelihood(
    response_values: list[float],
    residuals: list[float],
    inverse_covariance: list[list[float]],
    covariance: list[list[float]],
) -> float:
    sigma_squared = max(
        _quadratic_form(residuals, inverse_covariance) / len(response_values),
        1e-12,
    )
    return -0.5 * (
        len(response_values) * math.log(2.0 * math.pi * sigma_squared)
        + log_determinant(covariance)
        + len(response_values)
    )


def _estimate_lambda_for_pgls(
    *,
    log_likelihood_at_lambda: Callable[[float], float],
    null_log_likelihood: float,
    brownian_log_likelihood: float,
    profile_step: float = 0.01,
    confidence_interval_drop: float = 1.920729410347062,
    coarse_step: float = 0.05,
    fine_step: float = 0.005,
) -> PGLSLambdaFitReport:
    coarse_values = _grid_values(0.0, 1.0, coarse_step)
    coarse_best_lambda = max(coarse_values, key=log_likelihood_at_lambda)
    fine_values = _grid_values(
        max(0.0, coarse_best_lambda - coarse_step),
        min(1.0, coarse_best_lambda + coarse_step),
        fine_step,
    )
    fine_best_lambda = max(fine_values, key=log_likelihood_at_lambda)
    lower_bound = max(0.0, fine_best_lambda - fine_step)
    upper_bound = min(1.0, fine_best_lambda + fine_step)
    resolved_lambda = _maximize_bounded_log_likelihood(
        log_likelihood_at_lambda,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
    )
    best_log_likelihood = log_likelihood_at_lambda(resolved_lambda)
    threshold = best_log_likelihood - confidence_interval_drop
    profile_rows = [
        PGLSLambdaProfileRow(
            lambda_value=candidate,
            log_likelihood=log_likelihood_at_lambda(candidate),
            delta_log_likelihood=best_log_likelihood
            - log_likelihood_at_lambda(candidate),
            within_95_confidence_interval=(
                log_likelihood_at_lambda(candidate) >= threshold
            ),
        )
        for candidate in _grid_values(0.0, 1.0, profile_step)
    ]
    supported_rows = [row for row in profile_rows if row.within_95_confidence_interval]
    lower_bound = supported_rows[0].lambda_value if supported_rows else None
    upper_bound = supported_rows[-1].lambda_value if supported_rows else None
    return PGLSLambdaFitReport(
        mode="estimated",
        lambda_value=resolved_lambda,
        log_likelihood=best_log_likelihood,
        null_log_likelihood=null_log_likelihood,
        brownian_log_likelihood=brownian_log_likelihood,
        lower_95_confidence_interval=lower_bound,
        upper_95_confidence_interval=upper_bound,
        profile_rows=profile_rows,
    )


def _maximize_bounded_log_likelihood(
    log_likelihood_at_lambda: Callable[[float], float],
    *,
    lower_bound: float,
    upper_bound: float,
    tolerance: float = 1e-9,
    max_iterations: int = 80,
) -> float:
    if upper_bound <= lower_bound:
        return lower_bound
    inverse_phi = (math.sqrt(5.0) - 1.0) / 2.0
    left = lower_bound
    right = upper_bound
    interior_left = right - inverse_phi * (right - left)
    interior_right = left + inverse_phi * (right - left)
    left_value = log_likelihood_at_lambda(interior_left)
    right_value = log_likelihood_at_lambda(interior_right)
    for _ in range(max_iterations):
        if abs(right - left) <= tolerance:
            break
        if left_value <= right_value:
            left = interior_left
            interior_left = interior_right
            left_value = right_value
            interior_right = left + inverse_phi * (right - left)
            right_value = log_likelihood_at_lambda(interior_right)
        else:
            right = interior_right
            interior_right = interior_left
            right_value = left_value
            interior_left = right - inverse_phi * (right - left)
            left_value = log_likelihood_at_lambda(interior_left)
    candidate_values = {
        lower_bound: log_likelihood_at_lambda(lower_bound),
        upper_bound: log_likelihood_at_lambda(upper_bound),
        interior_left: left_value,
        interior_right: right_value,
        (left + right) / 2.0: log_likelihood_at_lambda((left + right) / 2.0),
    }
    return max(candidate_values, key=candidate_values.get)


def _lambda_log_likelihood(
    dataset: ComparativeDataset,
    design_matrix: list[list[float]],
    response_values: list[float],
    lambda_value: float,
) -> float:
    covariance = lambda_transform_covariance(dataset.covariance_matrix, lambda_value)
    inverse_covariance = invert_matrix(covariance)
    coefficients, _, fitted_values = _fit_gls(
        design_matrix,
        response_values,
        inverse_covariance,
    )
    del coefficients
    residuals = [
        observed - fitted
        for observed, fitted in zip(response_values, fitted_values, strict=True)
    ]
    return _gls_log_likelihood(
        response_values, residuals, inverse_covariance, covariance
    )


def _grid_values(start: float, stop: float, step: float) -> list[float]:
    values: list[float] = []
    current = start
    while current <= stop + (step / 2):
        values.append(round(current, 6))
        current += step
    return values
