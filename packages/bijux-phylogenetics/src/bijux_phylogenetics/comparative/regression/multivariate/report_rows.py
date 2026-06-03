from __future__ import annotations

from bijux_phylogenetics.comparative.pgls import PGLSResult

from .contracts import (
    MultivariateResponseCoefficientRow,
    MultivariateResponseModelRow,
)


def build_response_model_rows(
    response_models: list[PGLSResult],
) -> list[MultivariateResponseModelRow]:
    rows: list[MultivariateResponseModelRow] = []
    for model in response_models:
        residual_degrees_of_freedom = (
            model.coefficients[0].degrees_of_freedom if model.coefficients else 0
        )
        rows.append(
            MultivariateResponseModelRow(
                response=model.response,
                formula=model.formula.formula,
                predictor_term_count=len(model.predictors),
                encoded_term_count=len(model.encoded_columns),
                taxon_count=model.taxon_count,
                lambda_value=model.lambda_value,
                log_likelihood=model.log_likelihood,
                residual_variance=model.residual_variance,
                r_squared=model.r_squared,
                residual_degrees_of_freedom=residual_degrees_of_freedom,
            )
        )
    return rows


def build_response_coefficient_rows(
    response_models: list[PGLSResult],
) -> list[MultivariateResponseCoefficientRow]:
    rows: list[MultivariateResponseCoefficientRow] = []
    for model in response_models:
        for coefficient in model.coefficients:
            rows.append(
                MultivariateResponseCoefficientRow(
                    response=model.response,
                    formula=model.formula.formula,
                    term=coefficient.name,
                    estimate=coefficient.estimate,
                    standard_error=coefficient.standard_error,
                    test_statistic=coefficient.test_statistic,
                    p_value=coefficient.p_value,
                    lower_95_confidence_interval=coefficient.lower_95_confidence_interval,
                    upper_95_confidence_interval=coefficient.upper_95_confidence_interval,
                    degrees_of_freedom=coefficient.degrees_of_freedom,
                    inference_distribution=coefficient.inference_distribution,
                )
            )
    return rows


__all__ = [
    "build_response_coefficient_rows",
    "build_response_model_rows",
]
