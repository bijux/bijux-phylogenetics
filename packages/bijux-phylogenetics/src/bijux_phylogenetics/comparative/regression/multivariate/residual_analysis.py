from __future__ import annotations

import math

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    student_t_two_sided_p_value,
)
from bijux_phylogenetics.comparative.pgls import PGLSInputReport, PGLSResult
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .contracts import (
    MULTIVARIATE_LAMBDA_DIVERGENCE_WARNING_THRESHOLD,
    MULTIVARIATE_NEAR_SINGULAR_CONDITION_THRESHOLD,
    MULTIVARIATE_NUMERICAL_TOLERANCE,
    MULTIVARIATE_WEAK_SAMPLE_RESIDUAL_DF_THRESHOLD,
    MultivariateResidualAssociationRow,
    MultivariateResidualCorrelationRow,
    MultivariateResidualCovarianceDiagnostics,
    MultivariateResidualCovarianceRow,
)


def build_residual_covariance_rows(
    response_models: list[PGLSResult],
) -> list[MultivariateResidualCovarianceRow]:
    covariance_rows: list[MultivariateResidualCovarianceRow] = []
    for left_model in response_models:
        for right_model in response_models:
            covariance, correlation = covariance_and_correlation(
                left_model.residuals, right_model.residuals
            )
            covariance_rows.append(
                MultivariateResidualCovarianceRow(
                    left_response=left_model.response,
                    right_response=right_model.response,
                    covariance=covariance,
                    correlation=correlation,
                    pair_count=len(left_model.residuals),
                    is_diagonal=left_model.response == right_model.response,
                )
            )
    return covariance_rows


def build_residual_correlation_rows(
    covariance_rows: list[MultivariateResidualCovarianceRow],
) -> list[MultivariateResidualCorrelationRow]:
    return [
        MultivariateResidualCorrelationRow(
            left_response=row.left_response,
            right_response=row.right_response,
            correlation=row.correlation,
            pair_count=row.pair_count,
            is_diagonal=row.is_diagonal,
        )
        for row in covariance_rows
    ]


def build_residual_association_rows(
    response_models: list[PGLSResult],
) -> list[MultivariateResidualAssociationRow]:
    association_rows: list[MultivariateResidualAssociationRow] = []
    for left_index, left_model in enumerate(response_models):
        for right_model in response_models[left_index + 1 :]:
            covariance, correlation = covariance_and_correlation(
                left_model.residuals, right_model.residuals
            )
            pair_count = len(left_model.residuals)
            lower, upper = fisher_interval(correlation, pair_count)
            test_statistic, p_value = correlation_test(correlation, pair_count)
            association_rows.append(
                MultivariateResidualAssociationRow(
                    left_response=left_model.response,
                    right_response=right_model.response,
                    pair_count=pair_count,
                    covariance=covariance,
                    correlation=correlation,
                    test_statistic=test_statistic,
                    p_value=p_value,
                    lower_95_confidence_interval=lower,
                    upper_95_confidence_interval=upper,
                )
            )
    return association_rows


def build_residual_covariance_diagnostics(
    covariance_rows: list[MultivariateResidualCovarianceRow],
) -> MultivariateResidualCovarianceDiagnostics:
    response_names = ordered_response_names(covariance_rows)
    covariance_matrix = build_covariance_matrix(covariance_rows, response_names)
    response_count = len(response_names)
    rank = matrix_rank(covariance_matrix, tolerance=MULTIVARIATE_NUMERICAL_TOLERANCE)
    is_singular = rank < response_count
    if is_singular:
        condition_number = math.inf
    else:
        inverse = invert_matrix(covariance_matrix)
        condition_number = matrix_infinity_norm(
            covariance_matrix
        ) * matrix_infinity_norm(inverse)
    return MultivariateResidualCovarianceDiagnostics(
        response_count=response_count,
        matrix_rank=rank,
        condition_number=condition_number,
        is_singular=is_singular,
        is_near_singular=(
            is_singular
            or condition_number >= MULTIVARIATE_NEAR_SINGULAR_CONDITION_THRESHOLD
        ),
    )


def build_multivariate_warnings(
    *,
    final_input_reports: list[PGLSInputReport],
    response_models: list[PGLSResult],
    covariance_diagnostics: MultivariateResidualCovarianceDiagnostics,
) -> list[str]:
    warnings: list[str] = []
    for report in final_input_reports:
        for warning in report.warnings:
            warnings.append(f"{report.formula.response}: {warning}")
    minimum_residual_df = min(
        (
            model.coefficients[0].degrees_of_freedom
            for model in response_models
            if model.coefficients
        ),
        default=0,
    )
    if minimum_residual_df <= MULTIVARIATE_WEAK_SAMPLE_RESIDUAL_DF_THRESHOLD:
        warnings.append(
            "shared multivariate regression retains weak residual degrees of freedom, so response-specific coefficients and residual covariance estimates may be unstable"
        )
    lambda_span = response_lambda_span(response_models)
    if lambda_span >= MULTIVARIATE_LAMBDA_DIVERGENCE_WARNING_THRESHOLD:
        minimum_lambda = min(model.lambda_value for model in response_models)
        maximum_lambda = max(model.lambda_value for model in response_models)
        warnings.append(
            "response models resolved materially different Pagel lambda values "
            f"({format(minimum_lambda, '.15g')} to {format(maximum_lambda, '.15g')}), "
            "so shared residual covariance and correlation compare residuals fit under different phylogenetic error assumptions"
        )
    if covariance_diagnostics.is_near_singular:
        warnings.append(
            "residual covariance matrix is singular or near-singular within the multivariate numerical tolerance"
        )
    return deduplicate_preserving_order(warnings)


def deduplicate_preserving_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def response_lambda_span(response_models: list[PGLSResult]) -> float:
    if not response_models:
        return 0.0
    lambda_values = [model.lambda_value for model in response_models]
    return max(lambda_values) - min(lambda_values)


def covariance_and_correlation(
    left: list[float], right: list[float]
) -> tuple[float, float]:
    pair_count = len(left)
    if pair_count != len(right):
        raise ComparativeMethodError(
            "multivariate comparative regression requires aligned residual vectors"
        )
    if pair_count < 2:
        return 0.0, 0.0
    left_mean = sum(left) / pair_count
    right_mean = sum(right) / pair_count
    centered_pairs = [
        (left_value - left_mean, right_value - right_mean)
        for left_value, right_value in zip(left, right, strict=True)
    ]
    covariance = sum(x_value * y_value for x_value, y_value in centered_pairs) / (
        pair_count - 1
    )
    left_variance = sum(x_value * x_value for x_value, _ in centered_pairs) / (
        pair_count - 1
    )
    right_variance = sum(y_value * y_value for _, y_value in centered_pairs) / (
        pair_count - 1
    )
    denominator = math.sqrt(left_variance * right_variance)
    if denominator <= MULTIVARIATE_NUMERICAL_TOLERANCE:
        if all(
            math.isclose(
                left_value,
                right_value,
                rel_tol=0.0,
                abs_tol=MULTIVARIATE_NUMERICAL_TOLERANCE,
            )
            for left_value, right_value in zip(left, right, strict=True)
        ):
            return covariance, 1.0
        return covariance, 0.0
    correlation = covariance / denominator
    correlation = max(-1.0, min(1.0, correlation))
    return covariance, correlation


def ordered_response_names(
    covariance_rows: list[MultivariateResidualCovarianceRow],
) -> list[str]:
    response_names: list[str] = []
    for row in covariance_rows:
        if row.left_response not in response_names:
            response_names.append(row.left_response)
    return response_names


def build_covariance_matrix(
    covariance_rows: list[MultivariateResidualCovarianceRow],
    response_names: list[str],
) -> list[list[float]]:
    lookup = {
        (row.left_response, row.right_response): row.covariance
        for row in covariance_rows
    }
    return [
        [lookup[(left_response, right_response)] for right_response in response_names]
        for left_response in response_names
    ]


def matrix_rank(matrix: list[list[float]], *, tolerance: float) -> int:
    working = [list(row) for row in matrix]
    row_count = len(working)
    column_count = len(working[0]) if working else 0
    rank = 0
    pivot_row = 0
    for column_index in range(column_count):
        if pivot_row >= row_count:
            break
        best_row = max(
            range(pivot_row, row_count),
            key=lambda row_index: abs(working[row_index][column_index]),
        )
        pivot_value = working[best_row][column_index]
        if abs(pivot_value) <= tolerance:
            continue
        if best_row != pivot_row:
            working[pivot_row], working[best_row] = (
                working[best_row],
                working[pivot_row],
            )
        for row_index in range(pivot_row + 1, row_count):
            factor = working[row_index][column_index] / working[pivot_row][column_index]
            if abs(factor) <= tolerance:
                continue
            for trailing_index in range(column_index, column_count):
                working[row_index][trailing_index] -= (
                    factor * working[pivot_row][trailing_index]
                )
        rank += 1
        pivot_row += 1
    return rank


def matrix_infinity_norm(matrix: list[list[float]]) -> float:
    return max(
        (sum(abs(value) for value in row) for row in matrix),
        default=0.0,
    )


def correlation_test(correlation: float, pair_count: int) -> tuple[float, float]:
    if pair_count <= 2 or abs(correlation) >= 1.0:
        return math.inf if abs(correlation) >= 1.0 else 0.0, 0.0
    degrees_of_freedom = pair_count - 2
    test_statistic = correlation * math.sqrt(
        degrees_of_freedom
        / max(MULTIVARIATE_NUMERICAL_TOLERANCE, 1.0 - (correlation * correlation))
    )
    return test_statistic, student_t_two_sided_p_value(
        test_statistic, degrees_of_freedom
    )


def fisher_interval(
    correlation: float, pair_count: int
) -> tuple[float | None, float | None]:
    if pair_count <= 3 or abs(correlation) >= 1.0:
        return None, None
    fisher_z = math.atanh(correlation)
    standard_error = 1.0 / math.sqrt(pair_count - 3)
    lower = math.tanh(fisher_z - (1.959963984540054 * standard_error))
    upper = math.tanh(fisher_z + (1.959963984540054 * standard_error))
    return lower, upper


__all__ = [
    "build_multivariate_warnings",
    "build_residual_association_rows",
    "build_residual_correlation_rows",
    "build_residual_covariance_diagnostics",
    "build_residual_covariance_rows",
]
