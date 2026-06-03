from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    matrix_multiply,
    matrix_vector_multiply,
    stable_covariance,
    transpose,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    lambda_transform_covariance,
    load_comparative_dataset,
)
from bijux_phylogenetics.comparative.pgls import (
    ComparativeFormulaSpecification,
    PGLSFormulaAudit,
    inspect_pgls_inputs,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    MethodTierAssessment,
    phylogenetic_logistic_method_tier,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

_WALD_NORMAL_95_CRITICAL_VALUE = 1.959963984540054


@dataclass(slots=True)
class PhylogeneticLogisticWarning:
    """One reviewer-facing warning from the logistic approximation workflow."""

    code: str
    message: str


@dataclass(slots=True)
class PhylogeneticLogisticCoefficient:
    """One coefficient from a phylogenetic logistic approximation."""

    name: str
    estimate: float
    standard_error: float
    test_statistic: float
    p_value: float
    lower_95_confidence_interval: float
    upper_95_confidence_interval: float
    inference_distribution: str


@dataclass(slots=True)
class PhylogeneticLogisticFittedRow:
    """One observed-versus-fitted row from a phylogenetic logistic fit."""

    taxon: str
    observed_response: int
    fitted_probability: float
    linear_predictor: float
    residual: float


@dataclass(slots=True)
class PhylogeneticLogisticReport:
    """Reviewer-facing phylogenetic logistic approximation result."""

    tree_path: Path
    traits_path: Path
    response: str
    formula: ComparativeFormulaSpecification
    formula_audit: PGLSFormulaAudit
    approximation_method: str
    taxon_count: int
    success_count: int
    failure_count: int
    lambda_value: float
    converged: bool
    iteration_count: int
    binomial_log_likelihood: float
    separation_detected: bool
    method_tier: MethodTierAssessment
    warnings: list[PhylogeneticLogisticWarning]
    coefficients: list[PhylogeneticLogisticCoefficient]
    fitted_rows: list[PhylogeneticLogisticFittedRow]


def summarize_phylogenetic_logistic(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str | None = None,
    predictors: list[str] | None = None,
    formula: str | None = None,
    taxon_column: str | None = None,
    lambda_value: float = 1.0,
    max_iterations: int = 100,
    tolerance: float = 1e-8,
) -> PhylogeneticLogisticReport:
    """Fit a binary-response phylogenetic logistic approximation with a working correlation."""
    if not 0.0 <= lambda_value <= 1.0:
        raise ComparativeMethodError(
            "phylogenetic logistic working-correlation lambda must be between 0 and 1 inclusive"
        )

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

    taxa = list(input_report.analysis_taxa)
    response_values = [
        int(round(row.response_value)) for row in input_report.model_matrix.rows
    ]
    if any(value not in {0, 1} for value in response_values):
        raise ComparativeMethodError(
            "phylogenetic logistic regression requires a binary response encoded as 0 and 1"
        )
    success_count = sum(response_values)
    failure_count = len(response_values) - success_count
    if success_count == 0 or failure_count == 0:
        raise ComparativeMethodError(
            "phylogenetic logistic regression requires at least one success and one failure after pruning"
        )

    encoded_columns = list(input_report.model_matrix.encoded_columns)
    design_matrix = [
        [row.encoded_values[column] for column in encoded_columns]
        for row in input_report.model_matrix.rows
    ]
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=input_report.formula_audit.response_column,
        taxon_column=taxon_column,
        minimum_taxa=len(encoded_columns) + 1,
        require_rooted=True,
        require_binary=False,
    )
    dataset = ComparativeDataset(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        tree=dataset.tree,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        taxa=taxa,
        trait_values=[float(value) for value in response_values],
        covariance_matrix=_subset_covariance(
            dataset.covariance_matrix, dataset.taxa, taxa
        ),
        readiness=dataset.readiness,
    )
    correlation = _covariance_to_correlation(
        lambda_transform_covariance(dataset.covariance_matrix, lambda_value)
    )
    fit = _fit_phylogenetic_logistic(
        design_matrix=design_matrix,
        response_values=response_values,
        working_correlation=correlation,
        include_intercept=input_report.formula.include_intercept,
        max_iterations=max_iterations,
        tolerance=tolerance,
    )
    coefficients = _build_coefficient_rows(
        encoded_columns, fit.beta, fit.covariance_matrix
    )
    fitted_rows = [
        PhylogeneticLogisticFittedRow(
            taxon=taxon,
            observed_response=observed,
            fitted_probability=probability,
            linear_predictor=eta,
            residual=observed - probability,
        )
        for taxon, observed, probability, eta in zip(
            taxa,
            response_values,
            fit.fitted_probabilities,
            fit.linear_predictors,
            strict=True,
        )
    ]
    warnings = list(fit.warnings)
    if any(
        probability <= 1e-6 or probability >= 1.0 - 1e-6
        for probability in fit.fitted_probabilities
    ):
        warnings.append(
            PhylogeneticLogisticWarning(
                code="fitted_probability_boundary",
                message="one or more fitted probabilities are effectively on the 0/1 boundary, which is consistent with separation risk",
            )
        )
    if any(abs(value) >= 8.0 for value in fit.beta):
        warnings.append(
            PhylogeneticLogisticWarning(
                code="large_coefficient_magnitude",
                message="one or more fitted coefficients are large in magnitude, which is consistent with separation risk",
            )
        )
    warnings = _deduplicate_warnings(warnings)
    method_tier = phylogenetic_logistic_method_tier(
        "phylogenetic-working-correlation-gee"
    )
    return PhylogeneticLogisticReport(
        tree_path=tree_path,
        traits_path=traits_path,
        response=input_report.response,
        formula=input_report.formula,
        formula_audit=input_report.formula_audit,
        approximation_method="phylogenetic-working-correlation-gee",
        taxon_count=len(taxa),
        success_count=success_count,
        failure_count=failure_count,
        lambda_value=lambda_value,
        converged=fit.converged,
        iteration_count=fit.iteration_count,
        binomial_log_likelihood=_binomial_log_likelihood(
            response_values,
            fit.fitted_probabilities,
        ),
        separation_detected=bool(warnings),
        method_tier=method_tier,
        warnings=warnings,
        coefficients=coefficients,
        fitted_rows=fitted_rows,
    )


def write_phylogenetic_logistic_coefficient_table(
    path: Path, report: PhylogeneticLogisticReport
) -> Path:
    """Write one coefficient ledger for a phylogenetic logistic approximation."""
    return write_taxon_rows(
        path,
        columns=[
            "response",
            "term",
            "estimate",
            "standard_error",
            "test_statistic",
            "p_value",
            "lower_95_confidence_interval",
            "upper_95_confidence_interval",
            "inference_distribution",
            "approximation_method",
            "lambda_value",
            "taxon_count",
            "success_count",
            "failure_count",
            "converged",
            "iteration_count",
            "binomial_log_likelihood",
            "separation_detected",
        ],
        rows=[
            {
                "response": report.response,
                "term": row.name,
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
                "inference_distribution": row.inference_distribution,
                "approximation_method": report.approximation_method,
                "lambda_value": format(report.lambda_value, ".15g"),
                "taxon_count": str(report.taxon_count),
                "success_count": str(report.success_count),
                "failure_count": str(report.failure_count),
                "converged": str(report.converged).lower(),
                "iteration_count": str(report.iteration_count),
                "binomial_log_likelihood": format(
                    report.binomial_log_likelihood, ".15g"
                ),
                "separation_detected": str(report.separation_detected).lower(),
            }
            for row in report.coefficients
        ],
    )


def write_phylogenetic_logistic_fitted_table(
    path: Path, report: PhylogeneticLogisticReport
) -> Path:
    """Write one taxon-level fitted-probability ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "observed_response",
            "fitted_probability",
            "linear_predictor",
            "residual",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "observed_response": str(row.observed_response),
                "fitted_probability": format(row.fitted_probability, ".15g"),
                "linear_predictor": format(row.linear_predictor, ".15g"),
                "residual": format(row.residual, ".15g"),
            }
            for row in report.fitted_rows
        ],
    )


def write_phylogenetic_logistic_excluded_taxa_table(
    path: Path, report: PhylogeneticLogisticReport
) -> Path:
    """Write the explicit excluded-taxa ledger for one phylogenetic logistic request."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason", "details"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
                "details": row.details,
            }
            for row in report.formula_audit.excluded_taxa
        ],
    )


@dataclass(slots=True)
class _PhylogeneticLogisticFit:
    beta: list[float]
    covariance_matrix: list[list[float]]
    fitted_probabilities: list[float]
    linear_predictors: list[float]
    converged: bool
    iteration_count: int
    warnings: list[PhylogeneticLogisticWarning]


def _fit_phylogenetic_logistic(
    *,
    design_matrix: list[list[float]],
    response_values: list[int],
    working_correlation: list[list[float]],
    include_intercept: bool,
    max_iterations: int,
    tolerance: float,
) -> _PhylogeneticLogisticFit:
    beta = [0.0 for _ in range(len(design_matrix[0]))]
    if include_intercept:
        mean_response = sum(response_values) / len(response_values)
        clipped_mean = min(max(mean_response, 1e-6), 1.0 - 1e-6)
        beta[0] = math.log(clipped_mean / (1.0 - clipped_mean))
    correlation_inverse = invert_matrix(working_correlation)
    warnings: list[PhylogeneticLogisticWarning] = []
    converged = False
    iteration_count = 0
    for iteration in range(1, max_iterations + 1):
        iteration_count = iteration
        linear_predictors = _matrix_times_vector(design_matrix, beta)
        fitted_probabilities = [_sigmoid(value) for value in linear_predictors]
        variance_weights = [
            max(probability * (1.0 - probability), 1e-8)
            for probability in fitted_probabilities
        ]
        weighted_design = [
            [value * math.sqrt(weight) for value in row]
            for row, weight in zip(design_matrix, variance_weights, strict=True)
        ]
        scaled_residuals = [
            (observed - probability) / math.sqrt(weight)
            for observed, probability, weight in zip(
                response_values,
                fitted_probabilities,
                variance_weights,
                strict=True,
            )
        ]
        information_matrix = _xt_ax(weighted_design, correlation_inverse)
        try:
            information_inverse = invert_matrix(information_matrix)
        except ValueError:
            warnings.append(
                PhylogeneticLogisticWarning(
                    code="stabilized_information_matrix",
                    message="the working-information matrix was singular and required diagonal stabilization, which is consistent with separation or overfitting risk",
                )
            )
            information_matrix = stable_covariance(information_matrix, epsilon=1e-8)
            information_inverse = invert_matrix(information_matrix)
        score_vector = _xt_residual(
            weighted_design,
            correlation_inverse,
            scaled_residuals,
        )
        delta = _matrix_times_vector(information_inverse, score_vector)
        step_size = 1.0
        current_log_likelihood = _binomial_log_likelihood(
            response_values,
            fitted_probabilities,
        )
        while step_size >= 1e-6:
            candidate_beta = [
                coefficient + step_size * shift
                for coefficient, shift in zip(beta, delta, strict=True)
            ]
            candidate_linear_predictors = _matrix_times_vector(
                design_matrix, candidate_beta
            )
            candidate_probabilities = [
                _sigmoid(value) for value in candidate_linear_predictors
            ]
            candidate_log_likelihood = _binomial_log_likelihood(
                response_values,
                candidate_probabilities,
            )
            if candidate_log_likelihood >= current_log_likelihood - 1e-10:
                beta = candidate_beta
                linear_predictors = candidate_linear_predictors
                fitted_probabilities = candidate_probabilities
                break
            step_size /= 2.0
        max_delta = max(abs(step_size * value) for value in delta)
        if max_delta <= tolerance:
            converged = True
            covariance_matrix = information_inverse
            return _PhylogeneticLogisticFit(
                beta=beta,
                covariance_matrix=covariance_matrix,
                fitted_probabilities=fitted_probabilities,
                linear_predictors=linear_predictors,
                converged=True,
                iteration_count=iteration_count,
                warnings=warnings,
            )
    warnings.append(
        PhylogeneticLogisticWarning(
            code="convergence_limit_reached",
            message="the phylogenetic logistic approximation reached the iteration limit before the coefficient updates satisfied the convergence tolerance",
        )
    )
    final_linear_predictors = _matrix_times_vector(design_matrix, beta)
    final_probabilities = [_sigmoid(value) for value in final_linear_predictors]
    final_variance_weights = [
        max(probability * (1.0 - probability), 1e-8)
        for probability in final_probabilities
    ]
    final_weighted_design = [
        [value * math.sqrt(weight) for value in row]
        for row, weight in zip(design_matrix, final_variance_weights, strict=True)
    ]
    final_information = _xt_ax(final_weighted_design, correlation_inverse)
    final_information = stable_covariance(final_information, epsilon=1e-8)
    return _PhylogeneticLogisticFit(
        beta=beta,
        covariance_matrix=invert_matrix(final_information),
        fitted_probabilities=final_probabilities,
        linear_predictors=final_linear_predictors,
        converged=converged,
        iteration_count=iteration_count,
        warnings=warnings,
    )


def _build_coefficient_rows(
    encoded_columns: list[str],
    beta: list[float],
    covariance_matrix: list[list[float]],
) -> list[PhylogeneticLogisticCoefficient]:
    rows: list[PhylogeneticLogisticCoefficient] = []
    for index, name in enumerate(encoded_columns):
        standard_error = math.sqrt(max(covariance_matrix[index][index], 0.0))
        test_statistic = beta[index] / standard_error if standard_error else 0.0
        interval_radius = _WALD_NORMAL_95_CRITICAL_VALUE * standard_error
        rows.append(
            PhylogeneticLogisticCoefficient(
                name=name,
                estimate=beta[index],
                standard_error=standard_error,
                test_statistic=test_statistic,
                p_value=math.erfc(abs(test_statistic) / math.sqrt(2.0)),
                lower_95_confidence_interval=beta[index] - interval_radius,
                upper_95_confidence_interval=beta[index] + interval_radius,
                inference_distribution="wald-normal",
            )
        )
    return rows


def _covariance_to_correlation(
    covariance_matrix: list[list[float]],
) -> list[list[float]]:
    diagonal = [
        max(covariance_matrix[index][index], 1e-12)
        for index in range(len(covariance_matrix))
    ]
    return [
        [
            (
                1.0
                if row_index == column_index
                else covariance_matrix[row_index][column_index]
                / math.sqrt(diagonal[row_index] * diagonal[column_index])
            )
            for column_index in range(len(covariance_matrix))
        ]
        for row_index in range(len(covariance_matrix))
    ]


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


def _xt_ax(
    weighted_design: list[list[float]],
    correlation_inverse: list[list[float]],
) -> list[list[float]]:
    return matrix_multiply(
        transpose(weighted_design),
        matrix_multiply(correlation_inverse, weighted_design),
    )


def _xt_residual(
    weighted_design: list[list[float]],
    correlation_inverse: list[list[float]],
    scaled_residuals: list[float],
) -> list[float]:
    return _matrix_times_vector(
        transpose(weighted_design),
        _matrix_times_vector(correlation_inverse, scaled_residuals),
    )


def _matrix_times_vector(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return matrix_vector_multiply(matrix, vector)


def _sigmoid(value: float) -> float:
    if value >= 0.0:
        scale = math.exp(-value)
        return max(min(1.0 / (1.0 + scale), 1.0 - 1e-12), 1e-12)
    scale = math.exp(value)
    return max(min(scale / (1.0 + scale), 1.0 - 1e-12), 1e-12)


def _binomial_log_likelihood(
    response_values: list[int], fitted_probabilities: list[float]
) -> float:
    return sum(
        observed * math.log(probability) + (1 - observed) * math.log1p(-probability)
        for observed, probability in zip(
            response_values, fitted_probabilities, strict=True
        )
    )


def _deduplicate_warnings(
    warnings: list[PhylogeneticLogisticWarning],
) -> list[PhylogeneticLogisticWarning]:
    seen: set[tuple[str, str]] = set()
    deduplicated: list[PhylogeneticLogisticWarning] = []
    for warning in warnings:
        key = (warning.code, warning.message)
        if key in seen:
            continue
        seen.add(key)
        deduplicated.append(warning)
    return deduplicated
