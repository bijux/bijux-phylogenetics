from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import MultivariateComparativeRegressionReport


def write_multivariate_response_model_table(
    path: Path, report: MultivariateComparativeRegressionReport
) -> Path:
    """Write one response-model summary ledger for a multivariate fit."""
    return write_taxon_rows(
        path,
        columns=[
            "response",
            "formula",
            "predictor_term_count",
            "encoded_term_count",
            "taxon_count",
            "lambda_value",
            "log_likelihood",
            "residual_variance",
            "r_squared",
            "residual_degrees_of_freedom",
        ],
        rows=[
            {
                "response": row.response,
                "formula": row.formula,
                "predictor_term_count": row.predictor_term_count,
                "encoded_term_count": row.encoded_term_count,
                "taxon_count": row.taxon_count,
                "lambda_value": format(row.lambda_value, ".15g"),
                "log_likelihood": format(row.log_likelihood, ".15g"),
                "residual_variance": format(row.residual_variance, ".15g"),
                "r_squared": format(row.r_squared, ".15g"),
                "residual_degrees_of_freedom": row.residual_degrees_of_freedom,
            }
            for row in report.response_model_rows
        ],
    )


def write_multivariate_response_coefficient_table(
    path: Path, report: MultivariateComparativeRegressionReport
) -> Path:
    """Write one coefficient ledger across all responses in a multivariate fit."""
    return write_taxon_rows(
        path,
        columns=[
            "response",
            "formula",
            "term",
            "estimate",
            "standard_error",
            "test_statistic",
            "p_value",
            "lower_95_confidence_interval",
            "upper_95_confidence_interval",
            "degrees_of_freedom",
            "inference_distribution",
        ],
        rows=[
            {
                "response": row.response,
                "formula": row.formula,
                "term": row.term,
                "estimate": format(row.estimate, ".15g"),
                "standard_error": format(row.standard_error, ".15g"),
                "test_statistic": format(row.test_statistic, ".15g"),
                "p_value": format(row.p_value, ".15g"),
                "lower_95_confidence_interval": format(
                    row.lower_95_confidence_interval, ".15g"
                ),
                "upper_95_confidence_interval": format(
                    row.upper_95_confidence_interval, ".15g"
                ),
                "degrees_of_freedom": row.degrees_of_freedom,
                "inference_distribution": row.inference_distribution,
            }
            for row in report.coefficient_rows
        ],
    )


def write_multivariate_residual_covariance_table(
    path: Path, report: MultivariateComparativeRegressionReport
) -> Path:
    """Write one residual covariance ledger across response traits."""
    return write_taxon_rows(
        path,
        columns=[
            "left_response",
            "right_response",
            "pair_count",
            "is_diagonal",
            "covariance",
            "correlation",
        ],
        rows=[
            {
                "left_response": row.left_response,
                "right_response": row.right_response,
                "pair_count": row.pair_count,
                "is_diagonal": str(row.is_diagonal).lower(),
                "covariance": format(row.covariance, ".15g"),
                "correlation": format(row.correlation, ".15g"),
            }
            for row in report.covariance_rows
        ],
    )


def write_multivariate_residual_correlation_table(
    path: Path, report: MultivariateComparativeRegressionReport
) -> Path:
    """Write one residual correlation matrix ledger across response traits."""
    return write_taxon_rows(
        path,
        columns=[
            "left_response",
            "right_response",
            "pair_count",
            "is_diagonal",
            "correlation",
        ],
        rows=[
            {
                "left_response": row.left_response,
                "right_response": row.right_response,
                "pair_count": row.pair_count,
                "is_diagonal": str(row.is_diagonal).lower(),
                "correlation": format(row.correlation, ".15g"),
            }
            for row in report.correlation_rows
        ],
    )


def write_multivariate_residual_association_table(
    path: Path, report: MultivariateComparativeRegressionReport
) -> Path:
    """Write one residual trait-association ledger across response traits."""
    return write_taxon_rows(
        path,
        columns=[
            "left_response",
            "right_response",
            "pair_count",
            "covariance",
            "correlation",
            "test_statistic",
            "p_value",
            "lower_95_confidence_interval",
            "upper_95_confidence_interval",
        ],
        rows=[
            {
                "left_response": row.left_response,
                "right_response": row.right_response,
                "pair_count": row.pair_count,
                "covariance": format(row.covariance, ".15g"),
                "correlation": format(row.correlation, ".15g"),
                "test_statistic": format(row.test_statistic, ".15g"),
                "p_value": format(row.p_value, ".15g"),
                "lower_95_confidence_interval": (
                    ""
                    if row.lower_95_confidence_interval is None
                    else format(row.lower_95_confidence_interval, ".15g")
                ),
                "upper_95_confidence_interval": (
                    ""
                    if row.upper_95_confidence_interval is None
                    else format(row.upper_95_confidence_interval, ".15g")
                ),
            }
            for row in report.association_rows
        ],
    )


def write_multivariate_excluded_taxa_table(
    path: Path, report: MultivariateComparativeRegressionReport
) -> Path:
    """Write one explicit excluded-taxon ledger for multivariate fitting."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "reason",
            "missing_columns",
            "blocking_responses",
            "details",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
                "missing_columns": ",".join(row.missing_columns),
                "blocking_responses": ",".join(row.blocking_responses),
                "details": row.details,
            }
            for row in report.excluded_taxa
        ],
    )


__all__ = [
    "write_multivariate_excluded_taxa_table",
    "write_multivariate_residual_association_table",
    "write_multivariate_residual_correlation_table",
    "write_multivariate_residual_covariance_table",
    "write_multivariate_response_coefficient_table",
    "write_multivariate_response_model_table",
]
