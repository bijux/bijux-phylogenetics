from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class PGLSPredictorClassification:
    """Detected schema kind for one requested PGLS predictor."""

    name: str
    kind: str
    raw_term: str | None = None
    source_column: str | None = None
    transformation: str | None = None
    reference_level: str | None = None
    encoded_columns: list[str] | None = None
    observed_levels: list[str] | None = None
    level_counts: dict[str, int] | None = None


@dataclass(slots=True)
class ComparativeFormulaSpecification:
    """Auditable formula-style comparative model specification."""

    response: str
    formula: str
    predictors: list[str]
    interaction_terms: list[str]
    include_intercept: bool


@dataclass(slots=True)
class PGLSInteractionAudit:
    """Explicit expansion of one interaction term into encoded columns."""

    term: str
    component_terms: list[str]
    encoded_columns: list[str]


@dataclass(slots=True)
class PGLSTaxonExclusion:
    """One taxon excluded from comparative fitting and why."""

    taxon: str
    reason: str
    details: str


@dataclass(slots=True)
class PGLSFormulaAudit:
    """Reviewer-facing audit of the requested PGLS formula and exclusions."""

    response_term: str
    response_column: str
    predictor_terms: list[PGLSPredictorClassification]
    interaction_terms: list[PGLSInteractionAudit]
    transformed_terms: list[str]
    excluded_taxa: list[PGLSTaxonExclusion]
    includes_intercept: bool
    encoded_columns: list[str]
    analysis_taxa: list[str]
    parameter_count: int
    minimum_required_taxa: int
    residual_degrees_of_freedom: int
    overfit_guard_triggered: bool
    warnings: list[str]


@dataclass(slots=True)
class PGLSModelMatrixRow:
    """One taxon-level encoded row from a comparative formula design matrix."""

    taxon: str
    response_value: float
    encoded_values: dict[str, float]


@dataclass(slots=True)
class PGLSModelMatrixReport:
    """Reviewer-facing design matrix generated from one comparative formula."""

    formula: ComparativeFormulaSpecification
    response_column: str
    encoded_columns: list[str]
    row_count: int
    rows: list[PGLSModelMatrixRow]


@dataclass(slots=True)
class PGLSInputReport:
    """Method-specific readiness summary for a PGLS request."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    response: str
    formula: ComparativeFormulaSpecification
    predictors: list[PGLSPredictorClassification]
    formula_audit: PGLSFormulaAudit
    categorical_predictors: list[str]
    encoded_columns: list[str]
    analysis_taxa: list[str]
    residual_degrees_of_freedom: int
    model_matrix: PGLSModelMatrixReport
    ready: bool
    blockers: list[str]
    warnings: list[str]


@dataclass(slots=True)
class PGLSCoefficient:
    """One fitted PGLS regression coefficient."""

    name: str
    estimate: float
    standard_error: float
    test_statistic: float
    p_value: float
    lower_95_confidence_interval: float
    upper_95_confidence_interval: float
    degrees_of_freedom: int
    inference_distribution: str


@dataclass(slots=True)
class PGLSLambdaProfileRow:
    """One likelihood-profile row across candidate Pagel lambda values."""

    lambda_value: float
    log_likelihood: float
    delta_log_likelihood: float
    within_95_confidence_interval: bool


@dataclass(slots=True)
class PGLSLambdaFitReport:
    """Pagel lambda fit surface for one PGLS model."""

    mode: str
    lambda_value: float
    log_likelihood: float
    null_log_likelihood: float
    brownian_log_likelihood: float
    lower_95_confidence_interval: float | None
    upper_95_confidence_interval: float | None
    profile_rows: list[PGLSLambdaProfileRow]


@dataclass(slots=True)
class PGLSFittedObservation:
    """Observed-versus-fitted summary for one analyzed taxon."""

    taxon: str
    observed: float
    fitted: float
    residual: float


@dataclass(slots=True)
class PGLSLeverageRow:
    """Influence summary for one analyzed taxon."""

    taxon: str
    leverage: float
    standardized_residual: float


@dataclass(slots=True)
class PGLSResidualOutlier:
    """One taxon with a large standardized residual."""

    taxon: str
    residual: float
    standardized_residual: float


@dataclass(slots=True)
class PGLSDiagnosticsReport:
    """Residual and leverage diagnostics for a fitted PGLS model."""

    residual_mean: float
    leverage_rows: list[PGLSLeverageRow]
    outlier_taxa: list[PGLSResidualOutlier]
    fitted_observed_rows: list[PGLSFittedObservation]


@dataclass(slots=True)
class PGLSResult:
    """Generalized least-squares regression result over a phylogenetic covariance model."""

    tree_path: Path
    traits_path: Path
    response: str
    formula: ComparativeFormulaSpecification
    predictors: list[str]
    interaction_terms: list[str]
    encoded_columns: list[str]
    taxon_count: int
    lambda_value: float
    lambda_fit: PGLSLambdaFitReport
    log_likelihood: float
    aic: float
    residual_variance: float
    r_squared: float
    coefficients: list[PGLSCoefficient]
    fitted_values: list[float]
    residuals: list[float]
    taxa: list[str]
    diagnostics: PGLSDiagnosticsReport


@dataclass(slots=True)
class ComparativeHypothesisTestRow:
    """One coefficient-level significance test across many comparative fits."""

    response: str
    term: str
    estimate: float
    p_value: float
    adjusted_p_value: float
    significant: bool


@dataclass(slots=True)
class ComparativeMultipleTestingReport:
    """Multiple-testing correction across repeated PGLS analyses."""

    tree_path: Path
    traits_path: Path
    responses: list[str]
    predictors: list[str]
    adjustment_method: str
    family_size: int
    raw_significant_count: int
    adjusted_significant_count: int
    rows: list[ComparativeHypothesisTestRow]


@dataclass(frozen=True, slots=True)
class _FormulaTermDescriptor:
    raw_term: str
    source_column: str
    transformation: str | None
