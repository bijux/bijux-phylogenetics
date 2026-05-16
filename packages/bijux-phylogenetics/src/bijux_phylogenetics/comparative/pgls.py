from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
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
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.core.traits import validate_traits_table
from bijux_phylogenetics.runtime.errors import ComparativeMethodError
from bijux_phylogenetics.io.trees import load_tree


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


def inspect_pgls_inputs(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
) -> PGLSInputReport:
    """Inspect whether a PGLS request is valid for the given tree and trait table."""
    specification = _resolve_formula_specification(
        response=response,
        predictors=predictors,
        formula=formula,
    )
    response_descriptor = _parse_term_descriptor(specification.response)
    if response_descriptor.transformation is not None:
        raise ComparativeMethodError(
            "transformed response terms are not supported for PGLS"
        )
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=response_descriptor.source_column,
        taxon_column=taxon_column,
    )
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    tree = load_tree(tree_path)
    trait_report = validate_traits_table(traits_path, taxon_column=taxon_column)
    column_kinds = {column.name: column.kind for column in trait_report.trait_columns}
    blockers: list[str] = []
    warnings = list(readiness.warnings)

    if not readiness.rooted:
        blockers.append("PGLS requires a rooted tree")
    if not readiness.complete_branch_lengths:
        blockers.append("PGLS requires complete tree branch lengths")
    if response_descriptor.source_column not in column_kinds:
        blockers.append(
            f"trait table does not contain response column '{response_descriptor.source_column}'"
        )
    elif column_kinds[response_descriptor.source_column] != "numeric":
        blockers.append(
            f"response column '{response_descriptor.source_column}' must be numeric for PGLS"
        )

    predictor_reports: list[PGLSPredictorClassification] = []
    categorical_predictors: list[str] = []
    encoded_columns = ["intercept"] if specification.include_intercept else []
    for predictor in specification.predictors:
        term_descriptor = _parse_term_descriptor(predictor)
        kind = column_kinds.get(term_descriptor.source_column)
        if kind is None:
            blockers.append(
                f"trait table does not contain predictor column '{term_descriptor.source_column}'"
            )
            continue
        if kind != "numeric" and term_descriptor.transformation is not None:
            blockers.append(
                f"transformation '{term_descriptor.transformation}' requires numeric predictor column '{term_descriptor.source_column}'"
            )
            continue
        if kind == "numeric":
            if term_descriptor.transformation is not None:
                warnings.append(
                    f"predictor term '{predictor}' applies {term_descriptor.transformation} transformation to column '{term_descriptor.source_column}'"
                )
            predictor_reports.append(
                PGLSPredictorClassification(
                    name=predictor,
                    kind=kind,
                    raw_term=predictor,
                    source_column=term_descriptor.source_column,
                    transformation=term_descriptor.transformation,
                )
            )
            encoded_columns.append(predictor)
            continue

        categorical_predictors.append(predictor)
        levels = sorted(
            {
                row[term_descriptor.source_column]
                for row in table.rows
                if row[table.taxon_column] in tree.tip_names
                and row.get(term_descriptor.source_column, "")
            }
        )
        level_counts = {
            level: sum(
                1
                for row in table.rows
                if row[table.taxon_column] in tree.tip_names
                and row.get(term_descriptor.source_column, "") == level
            )
            for level in levels
        }
        if len(levels) < 2:
            blockers.append(
                f"categorical predictor '{term_descriptor.source_column}' requires at least two observed levels"
            )
            predictor_reports.append(
                PGLSPredictorClassification(
                    name=predictor,
                    kind=kind,
                    raw_term=predictor,
                    source_column=term_descriptor.source_column,
                    reference_level=levels[0] if levels else None,
                    encoded_columns=[],
                    observed_levels=levels,
                    level_counts=level_counts,
                )
            )
            continue
        reference_level = levels[0] if specification.include_intercept else None
        dummy_levels = levels[1:] if specification.include_intercept else levels
        dummy_columns = [f"{predictor}[{level}]" for level in dummy_levels]
        encoded_columns.extend(dummy_columns)
        if specification.include_intercept:
            warnings.append(
                f"categorical predictor '{predictor}' will be dummy-encoded with reference level '{reference_level}'"
            )
        else:
            warnings.append(
                f"categorical predictor '{predictor}' will be fully indicator-encoded because the formula excludes an intercept"
            )
        predictor_reports.append(
            PGLSPredictorClassification(
                name=predictor,
                kind=kind,
                raw_term=predictor,
                source_column=term_descriptor.source_column,
                reference_level=reference_level,
                encoded_columns=dummy_columns,
                observed_levels=levels,
                level_counts=level_counts,
            )
        )

    report_by_name = {report.name: report for report in predictor_reports}
    interaction_audits: list[PGLSInteractionAudit] = []
    for interaction in specification.interaction_terms:
        factor_names = interaction.split(":")
        missing_factors = [name for name in factor_names if name not in report_by_name]
        if missing_factors:
            blockers.append(
                f"interaction term '{interaction}' references unknown predictor(s): {', '.join(missing_factors)}"
            )
            continue
        interaction_columns = _interaction_column_names(interaction, report_by_name)
        encoded_columns.extend(interaction_columns)
        interaction_audits.append(
            PGLSInteractionAudit(
                term=interaction,
                component_terms=factor_names,
                encoded_columns=interaction_columns,
            )
        )

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    analysis_taxa: list[str] = []
    missing_tree_taxa: list[str] = []
    excluded_taxa: list[PGLSTaxonExclusion] = []
    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None:
            missing_tree_taxa.append(taxon)
            excluded_taxa.append(
                PGLSTaxonExclusion(
                    taxon=taxon,
                    reason="missing_from_trait_table",
                    details="taxon is present in the tree but absent from the trait table",
                )
            )
            continue
        required_columns = [response_descriptor.source_column] + [
            report.source_column or report.name for report in predictor_reports
        ]
        missing_columns = [
            column for column in required_columns if not row.get(column, "")
        ]
        if missing_columns:
            excluded_taxa.append(
                PGLSTaxonExclusion(
                    taxon=taxon,
                    reason="missing_value",
                    details=f"taxon is missing required value(s): {', '.join(sorted(set(missing_columns)))}",
                )
            )
            continue
        try:
            _coerce_numeric_value(
                row[response_descriptor.source_column], descriptor=response_descriptor
            )
            for predictor_report in predictor_reports:
                if predictor_report.kind == "numeric":
                    source_column = (
                        predictor_report.source_column or predictor_report.name
                    )
                    _coerce_numeric_value(
                        row[source_column],
                        descriptor=_parse_term_descriptor(
                            predictor_report.raw_term or predictor_report.name
                        ),
                    )
        except ValueError:
            excluded_taxa.append(
                PGLSTaxonExclusion(
                    taxon=taxon,
                    reason="non_numeric_or_invalid_value",
                    details="taxon has non-numeric or transformation-invalid value(s) required by the model",
                )
            )
            continue
        analysis_taxa.append(taxon)
    if missing_tree_taxa:
        blockers.append(
            "PGLS requires all analyzed taxa to be resolved against the trait table"
        )
    residual_degrees_of_freedom = len(analysis_taxa) - len(encoded_columns)
    if residual_degrees_of_freedom <= 0:
        blockers.append(
            "PGLS overfit guard requires at least one residual degree of freedom after predictor encoding"
        )
    elif residual_degrees_of_freedom == 1:
        warnings.append(
            "PGLS has only one residual degree of freedom after predictor encoding"
        )
    transformed_terms = sorted(
        [
            predictor_report.name
            for predictor_report in predictor_reports
            if predictor_report.transformation is not None
        ]
    )
    formula_audit = PGLSFormulaAudit(
        response_term=specification.response,
        response_column=response_descriptor.source_column,
        predictor_terms=predictor_reports,
        interaction_terms=interaction_audits,
        transformed_terms=transformed_terms,
        excluded_taxa=excluded_taxa,
        includes_intercept=specification.include_intercept,
        encoded_columns=encoded_columns,
        analysis_taxa=analysis_taxa,
        parameter_count=len(encoded_columns),
        minimum_required_taxa=len(encoded_columns) + 1,
        residual_degrees_of_freedom=residual_degrees_of_freedom,
        overfit_guard_triggered=residual_degrees_of_freedom <= 0,
        warnings=warnings,
    )
    model_matrix = _build_model_matrix_report(
        rows_by_taxon=rows_by_taxon,
        taxa=analysis_taxa,
        specification=specification,
        predictor_reports=predictor_reports,
        response_descriptor=response_descriptor,
        response_column=response_descriptor.source_column,
    )

    return PGLSInputReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        response=specification.response,
        formula=specification,
        predictors=predictor_reports,
        formula_audit=formula_audit,
        categorical_predictors=categorical_predictors,
        encoded_columns=encoded_columns,
        analysis_taxa=analysis_taxa,
        residual_degrees_of_freedom=residual_degrees_of_freedom,
        model_matrix=model_matrix,
        ready=not blockers,
        blockers=blockers,
        warnings=warnings,
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


def run_pgls_multiple_testing(
    tree_path: Path,
    traits_path: Path,
    *,
    responses: list[str],
    predictors: list[str],
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
    adjustment_method: str = "benjamini-hochberg",
) -> ComparativeMultipleTestingReport:
    """Run repeated PGLS fits and adjust coefficient p-values across the family of tests."""
    if not responses:
        raise ComparativeMethodError(
            "multiple-testing analysis requires at least one response trait"
        )
    if adjustment_method != "benjamini-hochberg":
        raise ComparativeMethodError(
            "supported multiple-testing adjustment is 'benjamini-hochberg'"
        )
    rows: list[ComparativeHypothesisTestRow] = []
    for response in responses:
        report = run_pgls(
            tree_path,
            traits_path,
            response=response,
            predictors=predictors,
            taxon_column=taxon_column,
            lambda_value=lambda_value,
        )
        for coefficient in report.coefficients:
            if coefficient.name == "intercept":
                continue
            rows.append(
                ComparativeHypothesisTestRow(
                    response=response,
                    term=coefficient.name,
                    estimate=coefficient.estimate,
                    p_value=coefficient.p_value,
                    adjusted_p_value=coefficient.p_value,
                    significant=False,
                )
            )
    adjusted = _benjamini_hochberg_adjustment([row.p_value for row in rows])
    for row, adjusted_p_value in zip(rows, adjusted, strict=True):
        row.adjusted_p_value = adjusted_p_value
        row.significant = adjusted_p_value <= 0.05
    raw_significant_count = sum(1 for row in rows if row.p_value <= 0.05)
    adjusted_significant_count = sum(1 for row in rows if row.significant)
    return ComparativeMultipleTestingReport(
        tree_path=tree_path,
        traits_path=traits_path,
        responses=list(responses),
        predictors=list(predictors),
        adjustment_method=adjustment_method,
        family_size=len(rows),
        raw_significant_count=raw_significant_count,
        adjusted_significant_count=adjusted_significant_count,
        rows=rows,
    )


def build_pgls_model_matrix(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
) -> PGLSModelMatrixReport:
    """Build the encoded design matrix implied by one PGLS request."""
    return inspect_pgls_inputs(
        tree_path,
        traits_path,
        response=response,
        predictors=predictors,
        formula=formula,
        taxon_column=taxon_column,
    ).model_matrix


def write_pgls_model_matrix_table(path: Path, report: PGLSModelMatrixReport) -> Path:
    """Write a comparative model matrix as CSV or TSV."""
    return write_taxon_rows(
        path,
        columns=["taxon", "response_value", *report.encoded_columns],
        rows=[
            {
                "taxon": row.taxon,
                "response_value": format(row.response_value, ".15g"),
                **{
                    column: format(row.encoded_values[column], ".15g")
                    for column in report.encoded_columns
                },
            }
            for row in report.rows
        ],
    )


def _build_design_matrix(
    rows_by_taxon: dict[str, dict[str, str]],
    taxa: list[str],
    predictors: list[str],
    predictor_reports: list[PGLSPredictorClassification],
    interaction_terms: list[str],
    *,
    include_intercept: bool,
) -> tuple[list[list[float]], list[str]]:
    encoded_columns = ["intercept"] if include_intercept else []
    report_by_name = {report.name: report for report in predictor_reports}
    for predictor in predictors:
        report = report_by_name[predictor]
        if report.kind == "numeric":
            encoded_columns.append(predictor)
        else:
            encoded_columns.extend(report.encoded_columns or [])
    interaction_columns = {
        interaction: _interaction_column_names(interaction, report_by_name)
        for interaction in interaction_terms
    }
    for columns in interaction_columns.values():
        encoded_columns.extend(columns)
    matrix: list[list[float]] = []
    for taxon in taxa:
        row = rows_by_taxon[taxon]
        encoded_row = [1.0] if include_intercept else []
        encoded_main_effects: dict[str, list[tuple[str, float]]] = {}
        for predictor in predictors:
            report = report_by_name[predictor]
            if report.kind == "numeric":
                source_column = report.source_column or predictor
                numeric_value = _coerce_numeric_value(
                    row[source_column],
                    descriptor=_parse_term_descriptor(report.raw_term or predictor),
                )
                encoded_row.append(numeric_value)
                encoded_main_effects[predictor] = [(predictor, numeric_value)]
                continue
            categorical_rows: list[tuple[str, float]] = []
            for encoded_name in report.encoded_columns or []:
                level = encoded_name.removeprefix(f"{predictor}[").removesuffix("]")
                source_column = report.source_column or predictor
                value = 1.0 if row[source_column] == level else 0.0
                categorical_rows.append((encoded_name, value))
                encoded_row.append(value)
            encoded_main_effects[predictor] = categorical_rows
        for interaction in interaction_terms:
            interaction_values = _interaction_values(
                interaction,
                encoded_main_effects,
            )
            encoded_row.extend(value for _, value in interaction_values)
        matrix.append(encoded_row)
    return matrix, encoded_columns


def _build_model_matrix_report(
    *,
    rows_by_taxon: dict[str, dict[str, str]],
    taxa: list[str],
    specification: ComparativeFormulaSpecification,
    predictor_reports: list[PGLSPredictorClassification],
    response_descriptor: _FormulaTermDescriptor,
    response_column: str,
) -> PGLSModelMatrixReport:
    design_matrix, encoded_columns = _build_design_matrix(
        rows_by_taxon,
        taxa,
        specification.predictors,
        predictor_reports,
        specification.interaction_terms,
        include_intercept=specification.include_intercept,
    )
    rows = [
        PGLSModelMatrixRow(
            taxon=taxon,
            response_value=_coerce_numeric_value(
                rows_by_taxon[taxon][response_column],
                descriptor=response_descriptor,
            ),
            encoded_values={
                column: design_matrix[row_index][column_index]
                for column_index, column in enumerate(encoded_columns)
            },
        )
        for row_index, taxon in enumerate(taxa)
    ]
    return PGLSModelMatrixReport(
        formula=specification,
        response_column=response_column,
        encoded_columns=encoded_columns,
        row_count=len(rows),
        rows=rows,
    )


def _parse_term_descriptor(raw_term: str) -> _FormulaTermDescriptor:
    raw_term = raw_term.strip()
    if "(" not in raw_term:
        return _FormulaTermDescriptor(
            raw_term=raw_term,
            source_column=raw_term,
            transformation=None,
        )
    if not raw_term.endswith(")") or raw_term.count("(") != 1:
        raise ComparativeMethodError(
            f"unsupported comparative term syntax '{raw_term}'"
        )
    transformation, inner = raw_term.split("(", 1)
    transformation = transformation.strip()
    source_column = inner[:-1].strip()
    if transformation not in {"log", "log10", "sqrt"}:
        raise ComparativeMethodError(
            f"unsupported comparative transformation '{transformation}' in term '{raw_term}'"
        )
    if not source_column:
        raise ComparativeMethodError(
            f"comparative transformation term '{raw_term}' is missing a source column"
        )
    return _FormulaTermDescriptor(
        raw_term=raw_term,
        source_column=source_column,
        transformation=transformation,
    )


def _parse_right_hand_side_terms(right_hand_side: str) -> list[tuple[str, str]]:
    terms: list[tuple[str, str]] = []
    current: list[str] = []
    current_sign = "+"
    parenthesis_depth = 0
    for character in right_hand_side:
        if character == "(":
            parenthesis_depth += 1
            current.append(character)
            continue
        if character == ")":
            parenthesis_depth -= 1
            if parenthesis_depth < 0:
                raise ComparativeMethodError(
                    "comparative formula has unmatched closing parenthesis"
                )
            current.append(character)
            continue
        if parenthesis_depth == 0 and character in {"+", "-"}:
            term = "".join(current).strip()
            if term:
                terms.append((current_sign, term))
            current = []
            current_sign = character
            continue
        current.append(character)
    if parenthesis_depth != 0:
        raise ComparativeMethodError(
            "comparative formula has unmatched opening parenthesis"
        )
    trailing_term = "".join(current).strip()
    if trailing_term:
        terms.append((current_sign, trailing_term))
    return terms


def _coerce_numeric_value(
    raw_value: str, *, descriptor: _FormulaTermDescriptor
) -> float:
    value = float(raw_value)
    if descriptor.transformation is None:
        return value
    if descriptor.transformation == "log":
        if value <= 0.0:
            raise ValueError("log transformation requires strictly positive values")
        return math.log(value)
    if descriptor.transformation == "log10":
        if value <= 0.0:
            raise ValueError("log10 transformation requires strictly positive values")
        return math.log10(value)
    if descriptor.transformation == "sqrt":
        if value < 0.0:
            raise ValueError("sqrt transformation requires non-negative values")
        return math.sqrt(value)
    raise ComparativeMethodError(
        f"unsupported transformation '{descriptor.transformation}'"
    )


def _resolve_formula_specification(
    *,
    response: str | None,
    predictors: list[str] | None,
    formula: str | None,
) -> ComparativeFormulaSpecification:
    if formula:
        if response is not None or predictors:
            raise ComparativeMethodError(
                "provide either a formula or explicit response/predictors, not both"
            )
        return _parse_formula(formula)
    if response is None:
        raise ComparativeMethodError(
            "PGLS requires a response column when no formula is provided"
        )
    requested_predictors = list(predictors or [])
    if not requested_predictors:
        raise ComparativeMethodError("PGLS requires at least one predictor column")
    return ComparativeFormulaSpecification(
        response=response,
        formula=f"{response} ~ {' + '.join(requested_predictors)}",
        predictors=requested_predictors,
        interaction_terms=[],
        include_intercept=True,
    )


def _parse_formula(formula: str) -> ComparativeFormulaSpecification:
    if "~" not in formula:
        raise ComparativeMethodError("comparative formula must contain '~'")
    response, right_hand_side = [part.strip() for part in formula.split("~", 1)]
    if not response:
        raise ComparativeMethodError(
            "comparative formula requires a response on the left-hand side"
        )
    raw_terms = _parse_right_hand_side_terms(right_hand_side)
    if not raw_terms:
        raise ComparativeMethodError(
            "comparative formula requires at least one predictor term"
        )
    predictors: list[str] = []
    interaction_terms: list[str] = []
    include_intercept = True
    for sign, raw_term in raw_terms:
        if raw_term in {"0", "1"}:
            if raw_term == "0" or sign == "-":
                include_intercept = False
            elif raw_term == "1":
                include_intercept = True
            continue
        if sign == "-":
            raise ComparativeMethodError(
                f"unsupported comparative formula subtraction for term '{raw_term}'"
            )
        if "*" in raw_term:
            factors = [
                factor.strip() for factor in raw_term.split("*") if factor.strip()
            ]
            if len(factors) < 2:
                raise ComparativeMethodError(
                    f"invalid interaction expansion '{raw_term}'"
                )
            for factor in factors:
                if factor not in predictors:
                    predictors.append(factor)
            interaction = ":".join(factors)
            if interaction not in interaction_terms:
                interaction_terms.append(interaction)
            continue
        if ":" in raw_term:
            factors = [
                factor.strip() for factor in raw_term.split(":") if factor.strip()
            ]
            if len(factors) < 2:
                raise ComparativeMethodError(f"invalid interaction term '{raw_term}'")
            interaction = ":".join(factors)
            if interaction not in interaction_terms:
                interaction_terms.append(interaction)
            continue
        if raw_term not in predictors:
            predictors.append(raw_term)
    if not predictors and not interaction_terms:
        raise ComparativeMethodError(
            "comparative formula requires at least one predictor term"
        )
    return ComparativeFormulaSpecification(
        response=response,
        formula=formula.strip(),
        predictors=predictors,
        interaction_terms=interaction_terms,
        include_intercept=include_intercept,
    )


def _interaction_column_names(
    interaction: str,
    predictor_reports: dict[str, PGLSPredictorClassification],
) -> list[str]:
    encoded_components: list[list[str]] = []
    for factor in interaction.split(":"):
        report = predictor_reports[factor]
        if report.kind == "numeric":
            encoded_components.append([factor])
        else:
            encoded_components.append(list(report.encoded_columns or []))
    column_names = [""]
    for component_names in encoded_components:
        column_names = [
            f"{left}:{right}".strip(":")
            for left in column_names
            for right in component_names
        ]
    return column_names


def _interaction_values(
    interaction: str,
    encoded_main_effects: dict[str, list[tuple[str, float]]],
) -> list[tuple[str, float]]:
    values = [("", 1.0)]
    for factor in interaction.split(":"):
        expanded: list[tuple[str, float]] = []
        for left_name, left_value in values:
            for right_name, right_value in encoded_main_effects[factor]:
                expanded.append(
                    (f"{left_name}:{right_name}".strip(":"), left_value * right_value)
                )
        values = expanded
    return values


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
    resolved_lambda = max(fine_values, key=log_likelihood_at_lambda)
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


def _benjamini_hochberg_adjustment(p_values: list[float]) -> list[float]:
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [1.0] * len(p_values)
    running_minimum = 1.0
    total = len(p_values)
    for rank, (index, p_value) in enumerate(reversed(indexed), start=1):
        denominator = total - rank + 1
        candidate = min(1.0, (p_value * total) / denominator)
        running_minimum = min(running_minimum, candidate)
        adjusted[index] = running_minimum
    return adjusted
