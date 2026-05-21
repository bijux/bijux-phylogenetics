from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative.pgls import (
    run_pgls,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .contracts import (
    MULTIVARIATE_MISSING_VALUE_POLICY,
    MULTIVARIATE_NEAR_SINGULAR_CONDITION_THRESHOLD,
    MULTIVARIATE_NUMERICAL_TOLERANCE,
    MULTIVARIATE_WEAK_SAMPLE_RESIDUAL_DF_THRESHOLD,
    MultivariateComparativeRegressionReport,
    MultivariateResidualAssociationRow,
    MultivariateResidualCorrelationRow,
    MultivariateResidualCovarianceDiagnostics,
    MultivariateResidualCovarianceRow,
    MultivariateResponseCoefficientRow,
    MultivariateResponseModelRow,
    MultivariateTaxonExclusion,
)
from .input_curation import (
    build_response_formula as _build_response_formula,
    build_shared_taxon_exclusions as _build_shared_taxon_exclusions,
    inspect_response_model as _inspect_response_model,
    raise_for_input_blockers as _raise_for_input_blockers,
    shared_analysis_taxa as _shared_analysis_taxa,
)
from .report_rows import build_response_coefficient_rows as _build_response_coefficient_rows
from .report_rows import build_response_model_rows as _build_response_model_rows
from .residual_analysis import (
    build_multivariate_warnings as _build_multivariate_warnings,
    build_residual_association_rows as _build_residual_association_rows,
    build_residual_correlation_rows as _build_residual_correlation_rows,
    build_residual_covariance_diagnostics as _build_residual_covariance_diagnostics,
    build_residual_covariance_rows as _build_residual_covariance_rows,
)


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

