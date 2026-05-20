from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    student_t_two_sided_p_value,
)
from bijux_phylogenetics.comparative.pgls import (
    PGLSInputReport,
    PGLSResult,
    PGLSTaxonExclusion,
    inspect_pgls_inputs,
    run_pgls,
)
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.runtime.errors import ComparativeMethodError
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree

MULTIVARIATE_NUMERICAL_TOLERANCE = 1e-12
MULTIVARIATE_MISSING_VALUE_POLICY = (
    "shared_complete_case_across_responses_and_predictor_terms"
)
MULTIVARIATE_WEAK_SAMPLE_RESIDUAL_DF_THRESHOLD = 2
MULTIVARIATE_NEAR_SINGULAR_CONDITION_THRESHOLD = 1e12


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


def run_multivariate_comparative_regression(
    tree_path: Path,
    traits_path: Path,
    *,
    responses: list[str],
    predictors: list[str],
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> MultivariateComparativeRegressionReport:
    """Fit one shared-taxon comparative workflow across multiple response traits."""
    if len(responses) < 2:
        raise ComparativeMethodError(
            "multivariate comparative regression requires at least two response traits"
        )
    if not predictors:
        raise ComparativeMethodError(
            "multivariate comparative regression requires at least one predictor"
        )

    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    tree = load_tree(tree_path)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    tree_taxa = set(tree.tip_names)
    table_taxa = set(table.taxa)
    excluded_taxa: list[MultivariateTaxonExclusion] = []
    required_terms = [*responses, *predictors]
    for taxon in sorted(tree_taxa - table_taxa):
        excluded_taxa.append(
            MultivariateTaxonExclusion(
                taxon=taxon,
                reason="missing_from_trait_table",
                missing_columns=list(required_terms),
                blocking_responses=list(responses),
                details="taxon is present in the tree but absent from the trait table",
            )
        )
    for taxon in sorted(table_taxa - tree_taxa):
        excluded_taxa.append(
            MultivariateTaxonExclusion(
                taxon=taxon,
                reason="missing_from_tree",
                missing_columns=[],
                blocking_responses=[],
                details="taxon is present in the trait table but absent from the tree",
            )
        )

    overlap_taxa = sorted(tree_taxa & table_taxa)
    if len(overlap_taxa) < 2:
        raise ComparativeMethodError(
            "multivariate comparative regression requires at least two taxa shared between the tree and trait table"
        )

    overlap_tree, _ = prune_tree_to_requested_taxa(tree_path, overlap_taxa)
    overlap_rows = [rows_by_taxon[taxon] for taxon in overlap_taxa]
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-multivariate-"
    ) as tmp_dir:
        overlap_tree_path = Path(tmp_dir) / "multivariate-overlap-tree.nwk"
        overlap_table_path = Path(tmp_dir) / "multivariate-overlap-traits.tsv"
        overlap_tree_path.write_text(
            dumps_newick(overlap_tree) + "\n", encoding="utf-8"
        )
        write_taxon_rows(
            overlap_table_path,
            columns=table.columns,
            rows=overlap_rows,
        )

        overlap_reports = [
            _inspect_response_model(
                overlap_tree_path,
                overlap_table_path,
                response=response,
                predictors=predictors,
                taxon_column=table.taxon_column,
            )
            for response in responses
        ]
        _raise_for_input_blockers(overlap_reports)

        shared_taxa = _shared_analysis_taxa(overlap_reports)
        excluded_taxa.extend(
            _build_shared_taxon_exclusions(
                overlap_taxa=overlap_taxa,
                responses=responses,
                overlap_reports=overlap_reports,
            )
        )
        if len(shared_taxa) < 2:
            raise ComparativeMethodError(
                "multivariate comparative regression does not retain enough shared complete-case taxa"
            )

        reduced_tree, _ = prune_tree_to_requested_taxa(overlap_tree_path, shared_taxa)
        reduced_rows = [rows_by_taxon[taxon] for taxon in shared_taxa]
        reduced_tree_path = Path(tmp_dir) / "multivariate-tree.nwk"
        reduced_table_path = Path(tmp_dir) / "multivariate-traits.tsv"
        reduced_tree_path.write_text(
            dumps_newick(reduced_tree) + "\n", encoding="utf-8"
        )
        write_taxon_rows(
            reduced_table_path,
            columns=table.columns,
            rows=reduced_rows,
        )

        final_input_reports = [
            _inspect_response_model(
                reduced_tree_path,
                reduced_table_path,
                response=response,
                predictors=predictors,
                taxon_column=table.taxon_column,
            )
            for response in responses
        ]
        _raise_for_input_blockers(final_input_reports)

        response_models = [
            run_pgls(
                reduced_tree_path,
                reduced_table_path,
                formula=_build_response_formula(response, predictors),
                taxon_column=table.taxon_column,
                lambda_value=lambda_value,
            )
            for response in responses
        ]

    for model in response_models:
        model.tree_path = tree_path
        model.traits_path = traits_path

    response_model_rows = _build_response_model_rows(response_models)
    coefficient_rows = _build_response_coefficient_rows(response_models)
    covariance_rows = _build_residual_covariance_rows(response_models)
    correlation_rows = _build_residual_correlation_rows(covariance_rows)
    association_rows = _build_residual_association_rows(response_models)
    covariance_diagnostics = _build_residual_covariance_diagnostics(covariance_rows)
    warnings = _build_multivariate_warnings(
        final_input_reports=final_input_reports,
        response_models=response_models,
        covariance_diagnostics=covariance_diagnostics,
    )
    return MultivariateComparativeRegressionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        responses=list(responses),
        predictors=list(predictors),
        taxon_column=table.taxon_column,
        missing_value_policy=MULTIVARIATE_MISSING_VALUE_POLICY,
        numerical_tolerance=MULTIVARIATE_NUMERICAL_TOLERANCE,
        analysis_taxa=list(response_models[0].taxa),
        excluded_taxa=excluded_taxa,
        response_models=response_models,
        response_model_rows=response_model_rows,
        coefficient_rows=coefficient_rows,
        covariance_rows=covariance_rows,
        correlation_rows=correlation_rows,
        association_rows=association_rows,
        covariance_diagnostics=covariance_diagnostics,
        warnings=warnings,
    )


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


def _build_response_formula(response: str, predictors: list[str]) -> str:
    return f"{response} ~ {' + '.join(predictors)}"


def _inspect_response_model(
    tree_path: Path,
    traits_path: Path,
    *,
    response: str,
    predictors: list[str],
    taxon_column: str,
) -> PGLSInputReport:
    return inspect_pgls_inputs(
        tree_path,
        traits_path,
        formula=_build_response_formula(response, predictors),
        taxon_column=taxon_column,
    )


def _raise_for_input_blockers(reports: list[PGLSInputReport]) -> None:
    blockers: list[str] = []
    for report in reports:
        for blocker in report.blockers:
            blockers.append(f"{report.formula.response}: {blocker}")
    if blockers:
        raise ComparativeMethodError("; ".join(blockers))


def _shared_analysis_taxa(reports: list[PGLSInputReport]) -> list[str]:
    if not reports:
        return []
    shared = set(reports[0].analysis_taxa)
    for report in reports[1:]:
        shared &= set(report.analysis_taxa)
    return sorted(shared)


def _build_shared_taxon_exclusions(
    *,
    overlap_taxa: list[str],
    responses: list[str],
    overlap_reports: list[PGLSInputReport],
) -> list[MultivariateTaxonExclusion]:
    exclusions_by_response: dict[str, dict[str, PGLSTaxonExclusion]] = {
        report.formula.response: {
            exclusion.taxon: exclusion
            for exclusion in report.formula_audit.excluded_taxa
        }
        for report in overlap_reports
    }
    excluded_rows: list[MultivariateTaxonExclusion] = []
    for taxon in overlap_taxa:
        missing_columns: set[str] = set()
        blocking_responses: list[str] = []
        invalid_details: list[str] = []
        missing_details: list[str] = []
        other_details: list[str] = []
        for response in responses:
            exclusion = exclusions_by_response.get(response, {}).get(taxon)
            if exclusion is None:
                continue
            blocking_responses.append(response)
            if exclusion.reason == "missing_value":
                parsed_missing = _parse_missing_columns(exclusion.details)
                missing_columns.update(parsed_missing)
                missing_details.append(exclusion.details)
            elif exclusion.reason == "non_numeric_or_invalid_value":
                invalid_details.append(exclusion.details)
            else:
                other_details.append(exclusion.details)
        if not blocking_responses:
            continue
        if invalid_details:
            reason = "invalid_required_values"
            details = "; ".join(sorted(set(invalid_details)))
        elif missing_columns:
            reason = "missing_required_values"
            details = "; ".join(sorted(set(missing_details)))
        else:
            reason = "excluded_from_shared_complete_case"
            details = "; ".join(sorted(set(other_details)))
        excluded_rows.append(
            MultivariateTaxonExclusion(
                taxon=taxon,
                reason=reason,
                missing_columns=sorted(missing_columns),
                blocking_responses=blocking_responses,
                details=details,
            )
        )
    return excluded_rows


def _parse_missing_columns(details: str) -> list[str]:
    prefix = "taxon is missing required value(s): "
    if not details.startswith(prefix):
        return []
    missing = details.removeprefix(prefix)
    return [column.strip() for column in missing.split(",") if column.strip()]


def _build_response_model_rows(
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


def _build_response_coefficient_rows(
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


def _build_residual_covariance_rows(
    response_models: list[PGLSResult],
) -> list[MultivariateResidualCovarianceRow]:
    covariance_rows: list[MultivariateResidualCovarianceRow] = []
    for left_model in response_models:
        for right_model in response_models:
            covariance, correlation = _covariance_and_correlation(
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


def _build_residual_correlation_rows(
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


def _build_residual_association_rows(
    response_models: list[PGLSResult],
) -> list[MultivariateResidualAssociationRow]:
    association_rows: list[MultivariateResidualAssociationRow] = []
    for left_index, left_model in enumerate(response_models):
        for right_model in response_models[left_index + 1 :]:
            covariance, correlation = _covariance_and_correlation(
                left_model.residuals, right_model.residuals
            )
            pair_count = len(left_model.residuals)
            lower, upper = _fisher_interval(correlation, pair_count)
            test_statistic, p_value = _correlation_test(correlation, pair_count)
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


def _build_residual_covariance_diagnostics(
    covariance_rows: list[MultivariateResidualCovarianceRow],
) -> MultivariateResidualCovarianceDiagnostics:
    response_names = _ordered_response_names(covariance_rows)
    covariance_matrix = _covariance_matrix(covariance_rows, response_names)
    response_count = len(response_names)
    matrix_rank = _matrix_rank(
        covariance_matrix,
        tolerance=MULTIVARIATE_NUMERICAL_TOLERANCE,
    )
    is_singular = matrix_rank < response_count
    if is_singular:
        condition_number = math.inf
    else:
        inverse = invert_matrix(covariance_matrix)
        condition_number = _matrix_infinity_norm(
            covariance_matrix
        ) * _matrix_infinity_norm(inverse)
    return MultivariateResidualCovarianceDiagnostics(
        response_count=response_count,
        matrix_rank=matrix_rank,
        condition_number=condition_number,
        is_singular=is_singular,
        is_near_singular=(
            is_singular
            or condition_number >= MULTIVARIATE_NEAR_SINGULAR_CONDITION_THRESHOLD
        ),
    )


def _build_multivariate_warnings(
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
    if covariance_diagnostics.is_near_singular:
        warnings.append(
            "residual covariance matrix is singular or near-singular within the multivariate numerical tolerance"
        )
    return _deduplicate_preserving_order(warnings)


def _deduplicate_preserving_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def _covariance_and_correlation(
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


def _ordered_response_names(
    covariance_rows: list[MultivariateResidualCovarianceRow],
) -> list[str]:
    response_names: list[str] = []
    for row in covariance_rows:
        if row.left_response not in response_names:
            response_names.append(row.left_response)
    return response_names


def _covariance_matrix(
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


def _matrix_rank(matrix: list[list[float]], *, tolerance: float) -> int:
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


def _matrix_infinity_norm(matrix: list[list[float]]) -> float:
    return max(
        (sum(abs(value) for value in row) for row in matrix),
        default=0.0,
    )


def _correlation_test(correlation: float, pair_count: int) -> tuple[float, float]:
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


def _fisher_interval(
    correlation: float, pair_count: int
) -> tuple[float | None, float | None]:
    if pair_count <= 3 or abs(correlation) >= 1.0:
        return None, None
    fisher_z = math.atanh(correlation)
    standard_error = 1.0 / math.sqrt(pair_count - 3)
    lower = math.tanh(fisher_z - (1.959963984540054 * standard_error))
    upper = math.tanh(fisher_z + (1.959963984540054 * standard_error))
    return lower, upper
