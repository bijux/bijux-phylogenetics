from __future__ import annotations

import math

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    quadratic_form,
)

from .contracts import CorrelatedTraitComparisonRow

_FISHER_95_Z = 1.959963984540054
_LOG_2PI = math.log(2.0 * math.pi)


def _estimate_trait_covariance(
    left_values: list[float],
    right_values: list[float],
) -> list[list[float]]:
    count = len(left_values)
    left_variance = sum(value * value for value in left_values) / count
    right_variance = sum(value * value for value in right_values) / count
    covariance = (
        sum(left * right for left, right in zip(left_values, right_values, strict=True))
        / count
    )
    return [
        [left_variance, covariance],
        [covariance, right_variance],
    ]


def _multivariate_normal_log_likelihood(
    *,
    observations: list[list[float]],
    covariance_matrix: list[list[float]],
) -> float:
    inverse_covariance = invert_matrix(covariance_matrix)
    log_det = log_determinant(covariance_matrix)
    dimension = len(covariance_matrix)
    total = 0.0
    for observation in observations:
        total += -0.5 * (
            (dimension * _LOG_2PI)
            + log_det
            + quadratic_form(observation, inverse_covariance)
        )
    return total


def _comparison_rows(
    *,
    independent_log_likelihood: float,
    independent_parameter_count: int,
    correlated_log_likelihood: float,
    correlated_parameter_count: int,
    independent_description: str,
    correlated_description: str,
) -> list[CorrelatedTraitComparisonRow]:
    independent_aic = _aic(independent_parameter_count, independent_log_likelihood)
    correlated_aic = _aic(correlated_parameter_count, correlated_log_likelihood)
    best_aic = min(independent_aic, correlated_aic)
    return [
        CorrelatedTraitComparisonRow(
            model_kind="independent",
            model_description=independent_description,
            parameter_count=independent_parameter_count,
            log_likelihood=independent_log_likelihood,
            aic=independent_aic,
            delta_aic=independent_aic - best_aic,
            selected=math.isclose(independent_aic, best_aic, abs_tol=1e-12),
        ),
        CorrelatedTraitComparisonRow(
            model_kind="correlated",
            model_description=correlated_description,
            parameter_count=correlated_parameter_count,
            log_likelihood=correlated_log_likelihood,
            aic=correlated_aic,
            delta_aic=correlated_aic - best_aic,
            selected=math.isclose(correlated_aic, best_aic, abs_tol=1e-12),
        ),
    ]


def _aic(parameter_count: int, log_likelihood: float) -> float:
    return (2.0 * parameter_count) - (2.0 * log_likelihood)


def _sample_covariance_and_correlation(
    left_values: list[float],
    right_values: list[float],
) -> tuple[float, float]:
    count = len(left_values)
    left_mean = sum(left_values) / count
    right_mean = sum(right_values) / count
    covariance = (
        sum(
            (left - left_mean) * (right - right_mean)
            for left, right in zip(left_values, right_values, strict=True)
        )
        / count
    )
    left_variance = sum((value - left_mean) ** 2 for value in left_values) / count
    right_variance = sum((value - right_mean) ** 2 for value in right_values) / count
    return covariance, _correlation(left_variance, right_variance, covariance)


def _correlation(
    left_variance: float,
    right_variance: float,
    covariance: float,
) -> float:
    denominator = math.sqrt(left_variance * right_variance)
    if math.isclose(denominator, 0.0, abs_tol=1e-12):
        return 0.0
    return covariance / denominator


def _fisher_interval(
    correlation: float,
    count: int,
) -> tuple[float | None, float | None]:
    if count <= 3 or abs(correlation) >= 1.0:
        return None, None
    fisher_z = math.atanh(correlation)
    standard_error = 1.0 / math.sqrt(count - 3)
    return (
        math.tanh(fisher_z - (_FISHER_95_Z * standard_error)),
        math.tanh(fisher_z + (_FISHER_95_Z * standard_error)),
    )


def _chi_square_survival(statistic: float, degrees_of_freedom: int) -> float:
    if statistic <= 0.0 or degrees_of_freedom <= 0:
        return 1.0
    if degrees_of_freedom == 1:
        return math.erfc(math.sqrt(statistic / 2.0))
    z_score = (
        ((statistic / degrees_of_freedom) ** (1.0 / 3.0))
        - (1.0 - (2.0 / (9.0 * degrees_of_freedom)))
    ) / math.sqrt(2.0 / (9.0 * degrees_of_freedom))
    return 0.5 * math.erfc(z_score / math.sqrt(2.0))


def _format_optional(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")
