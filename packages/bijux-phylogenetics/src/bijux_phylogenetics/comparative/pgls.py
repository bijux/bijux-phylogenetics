from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import invert_matrix, log_determinant, matrix_multiply, transpose
from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    lambda_transform_covariance,
    load_comparative_dataset,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.core.metadata import load_taxon_table
from bijux_phylogenetics.core.traits import validate_traits_table
from bijux_phylogenetics.errors import ComparativeMethodError


@dataclass(slots=True)
class PGLSPredictorClassification:
    """Detected schema kind for one requested PGLS predictor."""

    name: str
    kind: str
    reference_level: str | None = None
    encoded_columns: list[str] | None = None


@dataclass(slots=True)
class ComparativeFormulaSpecification:
    """Auditable formula-style comparative model specification."""

    response: str
    formula: str
    predictors: list[str]
    interaction_terms: list[str]


@dataclass(slots=True)
class PGLSInputReport:
    """Method-specific readiness summary for a PGLS request."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    response: str
    formula: ComparativeFormulaSpecification
    predictors: list[PGLSPredictorClassification]
    categorical_predictors: list[str]
    encoded_columns: list[str]
    analysis_taxa: list[str]
    residual_degrees_of_freedom: int
    ready: bool
    blockers: list[str]
    warnings: list[str]


@dataclass(slots=True)
class PGLSCoefficient:
    """One fitted PGLS regression coefficient."""

    name: str
    estimate: float
    standard_error: float
    z_score: float
    p_value: float


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
    log_likelihood: float
    residual_variance: float
    r_squared: float
    coefficients: list[PGLSCoefficient]
    fitted_values: list[float]
    residuals: list[float]
    taxa: list[str]
    diagnostics: PGLSDiagnosticsReport


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
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=specification.response,
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
    if specification.response not in column_kinds:
        blockers.append(f"trait table does not contain response column '{specification.response}'")
    elif column_kinds[specification.response] != "numeric":
        blockers.append(f"response column '{specification.response}' must be numeric for PGLS")

    predictor_reports: list[PGLSPredictorClassification] = []
    categorical_predictors: list[str] = []
    encoded_columns = ["intercept"]
    for predictor in specification.predictors:
        kind = column_kinds.get(predictor)
        if kind is None:
            blockers.append(f"trait table does not contain predictor column '{predictor}'")
            continue
        if kind == "numeric":
            predictor_reports.append(PGLSPredictorClassification(name=predictor, kind=kind))
            encoded_columns.append(predictor)
            continue

        categorical_predictors.append(predictor)
        levels = sorted(
            {
                row[predictor]
                for row in table.rows
                if row[table.taxon_column] in tree.tip_names and row.get(predictor, "")
            }
        )
        if len(levels) < 2:
            blockers.append(f"categorical predictor '{predictor}' requires at least two observed levels")
            predictor_reports.append(
                PGLSPredictorClassification(
                    name=predictor,
                    kind=kind,
                    reference_level=levels[0] if levels else None,
                    encoded_columns=[],
                )
            )
            continue
        reference_level = levels[0]
        dummy_columns = [f"{predictor}[{level}]" for level in levels[1:]]
        encoded_columns.extend(dummy_columns)
        warnings.append(
            f"categorical predictor '{predictor}' will be dummy-encoded with reference level '{reference_level}'"
        )
        predictor_reports.append(
            PGLSPredictorClassification(
                name=predictor,
                kind=kind,
                reference_level=reference_level,
                encoded_columns=dummy_columns,
                )
            )

    report_by_name = {report.name: report for report in predictor_reports}
    for interaction in specification.interaction_terms:
        factor_names = interaction.split(":")
        missing_factors = [name for name in factor_names if name not in report_by_name]
        if missing_factors:
            blockers.append(
                f"interaction term '{interaction}' references unknown predictor(s): {', '.join(missing_factors)}"
            )
            continue
        encoded_columns.extend(_interaction_column_names(interaction, report_by_name))

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    analysis_taxa: list[str] = []
    missing_tree_taxa: list[str] = []
    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None:
            missing_tree_taxa.append(taxon)
            continue
        requested_columns = [specification.response, *specification.predictors]
        if any(not row.get(column, "") for column in requested_columns):
            continue
        try:
            float(row[specification.response])
            for predictor_report in predictor_reports:
                if predictor_report.kind == "numeric":
                    float(row[predictor_report.name])
        except ValueError:
            continue
        analysis_taxa.append(taxon)
    if missing_tree_taxa:
        blockers.append("PGLS requires all analyzed taxa to be resolved against the trait table")
    residual_degrees_of_freedom = len(analysis_taxa) - len(encoded_columns)
    if residual_degrees_of_freedom <= 0:
        blockers.append(
            "PGLS overfit guard requires at least one residual degree of freedom after predictor encoding"
        )
    elif residual_degrees_of_freedom == 1:
        warnings.append("PGLS has only one residual degree of freedom after predictor encoding")

    return PGLSInputReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        response=specification.response,
        formula=specification,
        predictors=predictor_reports,
        categorical_predictors=categorical_predictors,
        encoded_columns=encoded_columns,
        analysis_taxa=analysis_taxa,
        residual_degrees_of_freedom=residual_degrees_of_freedom,
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
        trait=input_report.response,
        taxon_column=taxon_column,
        minimum_taxa=len(input_report.encoded_columns) + 1,
        require_rooted=True,
        require_binary=False,
    )
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    taxa = list(input_report.analysis_taxa)
    response_values = [float(rows_by_taxon[taxon][input_report.response]) for taxon in taxa]
    design_matrix, encoded_columns = _build_design_matrix(
        rows_by_taxon,
        taxa,
        input_report.formula.predictors,
        input_report.predictors,
        input_report.formula.interaction_terms,
    )
    dataset = ComparativeDataset(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        tree=dataset.tree,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        taxa=taxa,
        trait_values=response_values,
        covariance_matrix=_subset_covariance(dataset.covariance_matrix, dataset.taxa, taxa),
        readiness=dataset.readiness,
    )
    resolved_lambda = _resolve_lambda(dataset, lambda_value)
    covariance = lambda_transform_covariance(dataset.covariance_matrix, resolved_lambda)
    inverse_covariance = invert_matrix(covariance)
    coefficients, covariance_of_betas, fitted_values = _fit_gls(design_matrix, response_values, inverse_covariance)
    residuals = [observed - fitted for observed, fitted in zip(response_values, fitted_values, strict=True)]
    degrees_of_freedom = len(response_values) - len(coefficients)
    residual_variance = _quadratic_form(residuals, inverse_covariance) / degrees_of_freedom
    coefficient_reports: list[PGLSCoefficient] = []
    for index, name in enumerate(encoded_columns):
        standard_error = math.sqrt(max(covariance_of_betas[index][index] * residual_variance, 0.0))
        z_score = coefficients[index] / standard_error if standard_error else 0.0
        p_value = 2.0 * (1.0 - _normal_cdf(abs(z_score)))
        coefficient_reports.append(
            PGLSCoefficient(
                name=name,
                estimate=coefficients[index],
                standard_error=standard_error,
                z_score=z_score,
                p_value=p_value,
            )
        )
    mean_response = sum(response_values) / len(response_values)
    total_sum_of_squares = sum((value - mean_response) ** 2 for value in response_values)
    residual_sum_of_squares = sum(value * value for value in residuals)
    r_squared = 1.0 - (residual_sum_of_squares / total_sum_of_squares) if total_sum_of_squares else 1.0
    log_likelihood = _gls_log_likelihood(
        response_values,
        residuals,
        inverse_covariance,
        covariance,
    )
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
        log_likelihood=log_likelihood,
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


def _build_design_matrix(
    rows_by_taxon: dict[str, dict[str, str]],
    taxa: list[str],
    predictors: list[str],
    predictor_reports: list[PGLSPredictorClassification],
    interaction_terms: list[str],
) -> tuple[list[list[float]], list[str]]:
    encoded_columns = ["intercept"]
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
        encoded_row = [1.0]
        encoded_main_effects: dict[str, list[tuple[str, float]]] = {}
        for predictor in predictors:
            report = report_by_name[predictor]
            if report.kind == "numeric":
                numeric_value = float(row[predictor])
                encoded_row.append(numeric_value)
                encoded_main_effects[predictor] = [(predictor, numeric_value)]
                continue
            categorical_rows: list[tuple[str, float]] = []
            for encoded_name in report.encoded_columns or []:
                level = encoded_name.removeprefix(f"{predictor}[").removesuffix("]")
                value = 1.0 if row[predictor] == level else 0.0
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


def _resolve_formula_specification(
    *,
    response: str | None,
    predictors: list[str] | None,
    formula: str | None,
) -> ComparativeFormulaSpecification:
    if formula:
        if response is not None or predictors:
            raise ComparativeMethodError("provide either a formula or explicit response/predictors, not both")
        return _parse_formula(formula)
    if response is None:
        raise ComparativeMethodError("PGLS requires a response column when no formula is provided")
    requested_predictors = list(predictors or [])
    if not requested_predictors:
        raise ComparativeMethodError("PGLS requires at least one predictor column")
    return ComparativeFormulaSpecification(
        response=response,
        formula=f"{response} ~ {' + '.join(requested_predictors)}",
        predictors=requested_predictors,
        interaction_terms=[],
    )


def _parse_formula(formula: str) -> ComparativeFormulaSpecification:
    if "~" not in formula:
        raise ComparativeMethodError("comparative formula must contain '~'")
    response, right_hand_side = [part.strip() for part in formula.split("~", 1)]
    if not response:
        raise ComparativeMethodError("comparative formula requires a response on the left-hand side")
    raw_terms = [term.strip() for term in right_hand_side.split("+") if term.strip()]
    if not raw_terms:
        raise ComparativeMethodError("comparative formula requires at least one predictor term")
    predictors: list[str] = []
    interaction_terms: list[str] = []
    for raw_term in raw_terms:
        if "*" in raw_term:
            factors = [factor.strip() for factor in raw_term.split("*") if factor.strip()]
            if len(factors) < 2:
                raise ComparativeMethodError(f"invalid interaction expansion '{raw_term}'")
            for factor in factors:
                if factor not in predictors:
                    predictors.append(factor)
            interaction = ":".join(factors)
            if interaction not in interaction_terms:
                interaction_terms.append(interaction)
            continue
        if ":" in raw_term:
            factors = [factor.strip() for factor in raw_term.split(":") if factor.strip()]
            if len(factors) < 2:
                raise ComparativeMethodError(f"invalid interaction term '{raw_term}'")
            interaction = ":".join(factors)
            if interaction not in interaction_terms:
                interaction_terms.append(interaction)
            continue
        if raw_term not in predictors:
            predictors.append(raw_term)
    return ComparativeFormulaSpecification(
        response=response,
        formula=formula.strip(),
        predictors=predictors,
        interaction_terms=interaction_terms,
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
                expanded.append((f"{left_name}:{right_name}".strip(":"), left_value * right_value))
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
    xt_vinv_y = [sum(row[index] * response_values[index] for index in range(len(response_values))) for row in xt_vinv]
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


def _resolve_lambda(
    dataset: ComparativeDataset,
    lambda_value: float | str,
) -> float:
    if isinstance(lambda_value, (float, int)):
        if not 0.0 <= lambda_value <= 1.0:
            raise ComparativeMethodError("PGLS lambda must be between 0 and 1 inclusive")
        return float(lambda_value)
    if lambda_value != "estimate":
        raise ComparativeMethodError("PGLS lambda must be 'estimate' or a numeric value")
    return _estimate_lambda_for_dataset(dataset)


def _subset_covariance(
    covariance_matrix: list[list[float]],
    source_taxa: list[str],
    target_taxa: list[str],
) -> list[list[float]]:
    positions = {taxon: index for index, taxon in enumerate(source_taxa)}
    return [
        [covariance_matrix[positions[left_taxon]][positions[right_taxon]] for right_taxon in target_taxa]
        for left_taxon in target_taxa
    ]


def _quadratic_form(vector: list[float], matrix: list[list[float]]) -> float:
    total = 0.0
    for row_index, row in enumerate(matrix):
        total += vector[row_index] * sum(value * vector[column_index] for column_index, value in enumerate(row))
    return total


def _normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def _gls_log_likelihood(
    response_values: list[float],
    residuals: list[float],
    inverse_covariance: list[list[float]],
    covariance: list[list[float]],
) -> float:
    sigma_squared = _quadratic_form(residuals, inverse_covariance) / len(response_values)
    return -0.5 * (
        len(response_values) * math.log(2.0 * math.pi * sigma_squared)
        + log_determinant(covariance)
        + len(response_values)
    )


def _estimate_lambda_for_dataset(
    dataset: ComparativeDataset,
    *,
    coarse_step: float = 0.05,
    fine_step: float = 0.005,
) -> float:
    coarse_values = _grid_values(0.0, 1.0, coarse_step)
    coarse_best_lambda = max(coarse_values, key=lambda candidate: _lambda_log_likelihood(dataset, candidate))
    fine_values = _grid_values(
        max(0.0, coarse_best_lambda - coarse_step),
        min(1.0, coarse_best_lambda + coarse_step),
        fine_step,
    )
    return max(fine_values, key=lambda candidate: _lambda_log_likelihood(dataset, candidate))


def _lambda_log_likelihood(dataset: ComparativeDataset, lambda_value: float) -> float:
    covariance = lambda_transform_covariance(dataset.covariance_matrix, lambda_value)
    inverse_covariance = invert_matrix(covariance)
    coefficients, _, fitted_values = _fit_gls(
        [[1.0] for _ in dataset.taxa],
        dataset.trait_values,
        inverse_covariance,
    )
    del coefficients
    residuals = [
        observed - fitted
        for observed, fitted in zip(dataset.trait_values, fitted_values, strict=True)
    ]
    return _gls_log_likelihood(dataset.trait_values, residuals, inverse_covariance, covariance)


def _grid_values(start: float, stop: float, step: float) -> list[float]:
    values: list[float] = []
    current = start
    while current <= stop + (step / 2):
        values.append(round(current, 6))
        current += step
    return values
