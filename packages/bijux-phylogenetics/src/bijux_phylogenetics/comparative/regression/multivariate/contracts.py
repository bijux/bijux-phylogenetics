from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import PGLSResult

MULTIVARIATE_NUMERICAL_TOLERANCE = 1e-12
MULTIVARIATE_MISSING_VALUE_POLICY = (
    "shared_complete_case_across_responses_and_predictor_terms"
)
MULTIVARIATE_WEAK_SAMPLE_RESIDUAL_DF_THRESHOLD = 2
MULTIVARIATE_NEAR_SINGULAR_CONDITION_THRESHOLD = 1e12
MULTIVARIATE_LAMBDA_DIVERGENCE_WARNING_THRESHOLD = 0.2


@dataclass(slots=True)
class MultivariateTaxonExclusion:
    """One taxon excluded from a shared multivariate comparative analysis."""

    taxon: str
    reason: str
    missing_columns: list[str]
    blocking_responses: list[str]
    details: str


@dataclass(slots=True)
class MultivariateResidualCovarianceRow:
    """One pairwise residual covariance row between two responses."""

    left_response: str
    right_response: str
    covariance: float
    correlation: float
    pair_count: int
    is_diagonal: bool


@dataclass(slots=True)
class MultivariateResidualCorrelationRow:
    """One pairwise residual correlation row between two responses."""

    left_response: str
    right_response: str
    correlation: float
    pair_count: int
    is_diagonal: bool


@dataclass(slots=True)
class MultivariateResidualAssociationRow:
    """One pairwise residual-association test between two responses."""

    left_response: str
    right_response: str
    pair_count: int
    covariance: float
    correlation: float
    test_statistic: float
    p_value: float
    lower_95_confidence_interval: float | None
    upper_95_confidence_interval: float | None


@dataclass(slots=True)
class MultivariateResidualCovarianceDiagnostics:
    """Matrix-level diagnostics for multivariate residual covariance."""

    response_count: int
    matrix_rank: int
    condition_number: float
    is_singular: bool
    is_near_singular: bool


@dataclass(slots=True)
class MultivariateResponseCoefficientRow:
    """One coefficient row from one response model in a multivariate fit."""

    response: str
    formula: str
    term: str
    estimate: float
    standard_error: float
    test_statistic: float
    p_value: float
    lower_95_confidence_interval: float
    upper_95_confidence_interval: float
    degrees_of_freedom: int
    inference_distribution: str


@dataclass(slots=True)
class MultivariateResponseModelRow:
    """One response-level model summary from a multivariate fit."""

    response: str
    formula: str
    predictor_term_count: int
    encoded_term_count: int
    taxon_count: int
    lambda_value: float
    log_likelihood: float
    residual_variance: float
    r_squared: float
    residual_degrees_of_freedom: int


@dataclass(slots=True)
class MultivariateComparativeRegressionReport:
    """Shared-taxon multivariate comparative regression summary."""

    tree_path: Path
    traits_path: Path
    responses: list[str]
    predictors: list[str]
    taxon_column: str
    missing_value_policy: str
    numerical_tolerance: float
    analysis_taxa: list[str]
    excluded_taxa: list[MultivariateTaxonExclusion]
    response_models: list[PGLSResult]
    response_model_rows: list[MultivariateResponseModelRow]
    coefficient_rows: list[MultivariateResponseCoefficientRow]
    covariance_rows: list[MultivariateResidualCovarianceRow]
    correlation_rows: list[MultivariateResidualCorrelationRow]
    association_rows: list[MultivariateResidualAssociationRow]
    covariance_diagnostics: MultivariateResidualCovarianceDiagnostics
    warnings: list[str]


__all__ = [
    "MULTIVARIATE_MISSING_VALUE_POLICY",
    "MULTIVARIATE_LAMBDA_DIVERGENCE_WARNING_THRESHOLD",
    "MULTIVARIATE_NEAR_SINGULAR_CONDITION_THRESHOLD",
    "MULTIVARIATE_NUMERICAL_TOLERANCE",
    "MULTIVARIATE_WEAK_SAMPLE_RESIDUAL_DF_THRESHOLD",
    "MultivariateComparativeRegressionReport",
    "MultivariateResidualAssociationRow",
    "MultivariateResidualCorrelationRow",
    "MultivariateResidualCovarianceDiagnostics",
    "MultivariateResidualCovarianceRow",
    "MultivariateResponseCoefficientRow",
    "MultivariateResponseModelRow",
    "MultivariateTaxonExclusion",
]
