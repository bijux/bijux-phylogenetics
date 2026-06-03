from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.comparative.pgls import run_pgls
from bijux_phylogenetics.comparative.regression import (
    ComparativeRegressionModelExclusion,
    ComparativeRegressionModelRow,
    ComparativeRegressionModelSelectionReport,
    compare_comparative_regression_models,
    summarize_phylogenetic_logistic,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa

from .columns import comparative_model_table_columns
from .models import (
    SupplementaryComparativeModelRow,
    SupplementaryComparativeModelTableResult,
)
from .shared import stringify_list, write_dict_rows


def _serialize_comparative_model_row(
    row: SupplementaryComparativeModelRow,
) -> dict[str, object]:
    return {
        "tree_source": row.tree_source,
        "traits_source": row.traits_source,
        "response": row.response,
        "formula": row.formula,
        "model_family": row.model_family,
        "selected_criterion": row.selected_criterion,
        "best_formula": row.best_formula,
        "rank": row.rank,
        "selected": row.selected,
        "analysis_taxon_count": row.analysis_taxon_count,
        "excluded_taxon_count": row.excluded_taxon_count,
        "excluded_taxa": stringify_list(row.excluded_taxa),
        "encoded_columns": stringify_list(row.encoded_columns),
        "coefficient_name": row.coefficient_name,
        "estimate": row.estimate,
        "standard_error": row.standard_error,
        "test_statistic": row.test_statistic,
        "p_value": row.p_value,
        "lower_95_confidence_interval": row.lower_95_confidence_interval,
        "upper_95_confidence_interval": row.upper_95_confidence_interval,
        "inference_distribution": row.inference_distribution,
        "phylogenetic_parameter_name": row.phylogenetic_parameter_name,
        "phylogenetic_parameter_value": row.phylogenetic_parameter_value,
        "phylogenetic_parameter_estimated": row.phylogenetic_parameter_estimated,
        "log_likelihood": row.log_likelihood,
        "aic": row.aic,
        "aicc": row.aicc,
        "bic": row.bic,
        "delta_aicc": row.delta_aicc,
        "delta_bic": row.delta_bic,
        "akaike_weight": row.akaike_weight,
        "residual_mean": "" if row.residual_mean is None else row.residual_mean,
        "outlier_taxon_count": row.outlier_taxon_count,
        "outlier_taxa": stringify_list(row.outlier_taxa),
        "max_leverage": "" if row.max_leverage is None else row.max_leverage,
        "max_abs_standardized_residual": ""
        if row.max_abs_standardized_residual is None
        else row.max_abs_standardized_residual,
        "converged": "" if row.converged is None else row.converged,
        "iteration_count": "" if row.iteration_count is None else row.iteration_count,
        "separation_detected": row.separation_detected,
        "diagnostic_warning_count": row.diagnostic_warning_count,
        "diagnostic_warnings": stringify_list(row.diagnostic_warnings),
    }


def _write_comparative_model_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryComparativeModelRow],
) -> Path:
    return write_dict_rows(
        path,
        columns=columns,
        rows=[_serialize_comparative_model_row(row) for row in rows],
    )


def _serialize_comparative_exclusion(
    row: ComparativeRegressionModelExclusion,
) -> str:
    if row.missing_columns:
        return f"{row.taxon}:{row.reason}:{','.join(row.missing_columns)}"
    return f"{row.taxon}:{row.reason}"


def _prepare_shared_comparative_inputs(
    *,
    tree_path: Path,
    traits_path: Path,
    taxon_column: str | None,
    analysis_taxa: list[str],
) -> tuple[Path, Path]:
    reduced_tree, _ = prune_tree_to_requested_taxa(tree_path, analysis_taxa)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    reduced_rows = [rows_by_taxon[taxon] for taxon in analysis_taxa]
    temp_dir = Path(tempfile.mkdtemp(prefix="bijux-comparative-supplement-"))
    reduced_tree_path = temp_dir / "comparative-model-tree.nwk"
    reduced_table_path = temp_dir / "comparative-model-traits.tsv"
    reduced_tree_path.write_text(dumps_newick(reduced_tree) + "\n", encoding="utf-8")
    write_taxon_rows(reduced_table_path, columns=table.columns, rows=reduced_rows)
    return reduced_tree_path, reduced_table_path


def _build_pgls_comparative_rows(
    *,
    report: ComparativeRegressionModelSelectionReport,
    model_row: ComparativeRegressionModelRow,
    reduced_tree_path: Path,
    reduced_table_path: Path,
    taxon_column: str,
    lambda_value: float | str,
) -> list[SupplementaryComparativeModelRow]:
    fitted = run_pgls(
        reduced_tree_path,
        reduced_table_path,
        formula=model_row.formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    outlier_taxa = [row.taxon for row in fitted.diagnostics.outlier_taxa]
    max_leverage = max(
        (row.leverage for row in fitted.diagnostics.leverage_rows),
        default=None,
    )
    max_abs_standardized_residual = max(
        (abs(row.standardized_residual) for row in fitted.diagnostics.leverage_rows),
        default=None,
    )
    excluded_taxa = [
        _serialize_comparative_exclusion(row) for row in report.excluded_taxa
    ]
    return [
        SupplementaryComparativeModelRow(
            tree_source=str(report.tree_path),
            traits_source=str(report.traits_path),
            response=report.response,
            formula=model_row.formula,
            model_family=model_row.model_family,
            selected_criterion=report.selected_criterion,
            best_formula=report.best_formula,
            rank=model_row.rank,
            selected=model_row.selected,
            analysis_taxon_count=len(report.analysis_taxa),
            excluded_taxon_count=len(report.excluded_taxa),
            excluded_taxa=excluded_taxa,
            encoded_columns=list(model_row.encoded_columns),
            coefficient_name=coefficient.name,
            estimate=coefficient.estimate,
            standard_error=coefficient.standard_error,
            test_statistic=coefficient.test_statistic,
            p_value=coefficient.p_value,
            lower_95_confidence_interval=coefficient.lower_95_confidence_interval,
            upper_95_confidence_interval=coefficient.upper_95_confidence_interval,
            inference_distribution=coefficient.inference_distribution,
            phylogenetic_parameter_name=model_row.phylogenetic_parameter_name,
            phylogenetic_parameter_value=fitted.lambda_fit.lambda_value,
            phylogenetic_parameter_estimated=fitted.lambda_fit.mode == "estimated",
            log_likelihood=fitted.log_likelihood,
            aic=model_row.aic,
            aicc=model_row.aicc,
            bic=model_row.bic,
            delta_aicc=model_row.delta_aicc,
            delta_bic=model_row.delta_bic,
            akaike_weight=model_row.akaike_weight,
            residual_mean=fitted.diagnostics.residual_mean,
            outlier_taxon_count=len(outlier_taxa),
            outlier_taxa=outlier_taxa,
            max_leverage=max_leverage,
            max_abs_standardized_residual=max_abs_standardized_residual,
            converged=None,
            iteration_count=None,
            separation_detected=False,
            diagnostic_warning_count=0,
            diagnostic_warnings=[],
        )
        for coefficient in fitted.coefficients
    ]


def _build_logistic_comparative_rows(
    *,
    report: ComparativeRegressionModelSelectionReport,
    model_row: ComparativeRegressionModelRow,
    reduced_tree_path: Path,
    reduced_table_path: Path,
    taxon_column: str,
    lambda_value: float,
) -> list[SupplementaryComparativeModelRow]:
    fitted = summarize_phylogenetic_logistic(
        reduced_tree_path,
        reduced_table_path,
        formula=model_row.formula,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    residual_mean = (
        None
        if not fitted.fitted_rows
        else sum(row.residual for row in fitted.fitted_rows) / len(fitted.fitted_rows)
    )
    excluded_taxa = [
        _serialize_comparative_exclusion(row) for row in report.excluded_taxa
    ]
    warning_messages = [warning.message for warning in fitted.warnings]
    return [
        SupplementaryComparativeModelRow(
            tree_source=str(report.tree_path),
            traits_source=str(report.traits_path),
            response=report.response,
            formula=model_row.formula,
            model_family=model_row.model_family,
            selected_criterion=report.selected_criterion,
            best_formula=report.best_formula,
            rank=model_row.rank,
            selected=model_row.selected,
            analysis_taxon_count=len(report.analysis_taxa),
            excluded_taxon_count=len(report.excluded_taxa),
            excluded_taxa=excluded_taxa,
            encoded_columns=list(model_row.encoded_columns),
            coefficient_name=coefficient.name,
            estimate=coefficient.estimate,
            standard_error=coefficient.standard_error,
            test_statistic=coefficient.test_statistic,
            p_value=coefficient.p_value,
            lower_95_confidence_interval=coefficient.lower_95_confidence_interval,
            upper_95_confidence_interval=coefficient.upper_95_confidence_interval,
            inference_distribution=coefficient.inference_distribution,
            phylogenetic_parameter_name=model_row.phylogenetic_parameter_name,
            phylogenetic_parameter_value=fitted.lambda_value,
            phylogenetic_parameter_estimated=False,
            log_likelihood=fitted.binomial_log_likelihood,
            aic=model_row.aic,
            aicc=model_row.aicc,
            bic=model_row.bic,
            delta_aicc=model_row.delta_aicc,
            delta_bic=model_row.delta_bic,
            akaike_weight=model_row.akaike_weight,
            residual_mean=residual_mean,
            outlier_taxon_count=0,
            outlier_taxa=[],
            max_leverage=None,
            max_abs_standardized_residual=None,
            converged=fitted.converged,
            iteration_count=fitted.iteration_count,
            separation_detected=fitted.separation_detected,
            diagnostic_warning_count=len(fitted.warnings),
            diagnostic_warnings=warning_messages,
        )
        for coefficient in fitted.coefficients
    ]


def write_supplementary_comparative_model_table(
    path: Path,
    *,
    tree_path: Path,
    traits_path: Path,
    formulas: list[str],
    taxon_column: str | None = None,
    lambda_value: float | str = "estimate",
) -> SupplementaryComparativeModelTableResult:
    """Write one coefficient-level supplementary comparative-model table."""
    report = compare_comparative_regression_models(
        tree_path,
        traits_path,
        formulas=formulas,
        taxon_column=taxon_column,
        lambda_value=lambda_value,
    )
    reduced_tree_path, reduced_table_path = _prepare_shared_comparative_inputs(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=taxon_column,
        analysis_taxa=report.analysis_taxa,
    )
    rows: list[SupplementaryComparativeModelRow] = []
    for model_row in report.rows:
        if model_row.model_family == "pgls":
            rows.extend(
                _build_pgls_comparative_rows(
                    report=report,
                    model_row=model_row,
                    reduced_tree_path=reduced_tree_path,
                    reduced_table_path=reduced_table_path,
                    taxon_column=report.taxon_column,
                    lambda_value=lambda_value,
                )
            )
            continue
        rows.extend(
            _build_logistic_comparative_rows(
                report=report,
                model_row=model_row,
                reduced_tree_path=reduced_tree_path,
                reduced_table_path=reduced_table_path,
                taxon_column=report.taxon_column,
                lambda_value=float(lambda_value),
            )
        )
    columns = comparative_model_table_columns()
    _write_comparative_model_rows(path, columns=columns, rows=rows)
    return SupplementaryComparativeModelTableResult(
        output_path=path,
        row_count=len(rows),
        model_count=len(report.rows),
        selected_formula=report.best_formula,
        selected_criterion=report.selected_criterion,
        excluded_taxon_count=len(report.excluded_taxa),
        columns=columns,
        rows=rows,
    )
