from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative._math import student_t_two_sided_p_value
from bijux_phylogenetics.comparative.pgls import PGLSResult, run_pgls
from bijux_phylogenetics.core.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.errors import ComparativeMethodError
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class MultivariateTaxonExclusion:
    """One taxon excluded from a shared multivariate comparative analysis."""

    taxon: str
    reason: str
    missing_columns: list[str]


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
class MultivariateComparativeRegressionReport:
    """Shared-taxon multivariate comparative regression summary."""

    tree_path: Path
    traits_path: Path
    responses: list[str]
    predictors: list[str]
    taxon_column: str
    analysis_taxa: list[str]
    excluded_taxa: list[MultivariateTaxonExclusion]
    response_models: list[PGLSResult]
    covariance_rows: list[MultivariateResidualCovarianceRow]
    association_rows: list[MultivariateResidualAssociationRow]


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
    required_columns = [*responses, *predictors]
    for column in required_columns:
        if column not in table.columns:
            raise ComparativeMethodError(
                f"trait table does not contain required column '{column}'"
            )

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    tree_taxa = set(tree.tip_names)
    table_taxa = set(table.taxa)
    excluded_taxa: list[MultivariateTaxonExclusion] = []
    for taxon in sorted(tree_taxa - table_taxa):
        excluded_taxa.append(
            MultivariateTaxonExclusion(
                taxon=taxon,
                reason="missing_from_trait_table",
                missing_columns=list(required_columns),
            )
        )
    for taxon in sorted(table_taxa - tree_taxa):
        excluded_taxa.append(
            MultivariateTaxonExclusion(
                taxon=taxon,
                reason="missing_from_tree",
                missing_columns=[],
            )
        )

    overlap_taxa = sorted(tree_taxa & table_taxa)
    analysis_taxa: list[str] = []
    for taxon in overlap_taxa:
        row = rows_by_taxon[taxon]
        missing_columns = [column for column in required_columns if not row[column]]
        if missing_columns:
            excluded_taxa.append(
                MultivariateTaxonExclusion(
                    taxon=taxon,
                    reason="missing_required_values",
                    missing_columns=missing_columns,
                )
            )
            continue
        analysis_taxa.append(taxon)

    if len(analysis_taxa) < len(predictors) + 3:
        raise ComparativeMethodError(
            "multivariate comparative regression does not retain enough complete-case taxa"
        )

    reduced_tree, _ = prune_tree_to_requested_taxa(tree_path, analysis_taxa)
    reduced_rows = [rows_by_taxon[taxon] for taxon in analysis_taxa]
    with tempfile.TemporaryDirectory(prefix="bijux-phylogenetics-multivariate-") as tmp_dir:
        reduced_tree_path = Path(tmp_dir) / "multivariate-tree.nwk"
        reduced_table_path = Path(tmp_dir) / "multivariate-traits.tsv"
        reduced_tree_path.write_text(
            dumps_newick(reduced_tree) + "\n",
            encoding="utf-8",
        )
        write_taxon_rows(
            reduced_table_path,
            columns=table.columns,
            rows=reduced_rows,
        )
        response_models = [
            run_pgls(
                reduced_tree_path,
                reduced_table_path,
                response=response,
                predictors=predictors,
                taxon_column=table.taxon_column,
                lambda_value=lambda_value,
            )
            for response in responses
        ]
    for model in response_models:
        model.tree_path = tree_path
        model.traits_path = traits_path

    covariance_rows = _build_residual_covariance_rows(response_models)
    association_rows = _build_residual_association_rows(response_models)
    return MultivariateComparativeRegressionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        responses=list(responses),
        predictors=list(predictors),
        taxon_column=table.taxon_column,
        analysis_taxa=list(response_models[0].taxa),
        excluded_taxa=excluded_taxa,
        response_models=response_models,
        covariance_rows=covariance_rows,
        association_rows=association_rows,
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
        columns=["taxon", "reason", "missing_columns"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
                "missing_columns": ",".join(row.missing_columns),
            }
            for row in report.excluded_taxa
        ],
    )


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


def _covariance_and_correlation(
    left: list[float], right: list[float]
) -> tuple[float, float]:
    pair_count = len(left)
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
    correlation = covariance / denominator if denominator else 0.0
    return covariance, correlation


def _correlation_test(correlation: float, pair_count: int) -> tuple[float, float]:
    if pair_count <= 2 or abs(correlation) >= 1.0:
        return math.inf if abs(correlation) >= 1.0 else 0.0, 0.0
    degrees_of_freedom = pair_count - 2
    test_statistic = correlation * math.sqrt(
        degrees_of_freedom / max(1e-12, 1.0 - (correlation * correlation))
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
