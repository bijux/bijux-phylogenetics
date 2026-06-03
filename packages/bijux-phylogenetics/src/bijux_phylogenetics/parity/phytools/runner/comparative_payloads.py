from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.pgls import (
    build_pgls_model_matrix,
    inspect_pgls_inputs,
    run_pgls,
)
from bijux_phylogenetics.comparative.regression import (
    summarize_phylogenetic_anova,
    summarize_phylogenetic_residuals,
)

from ..registry import PhytoolsParityCase


def _pgls_parity_rows(*, model_matrix, report) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for coefficient in report.coefficients:
        rows.extend(
            [
                {
                    "row_kind": "coefficient_estimate",
                    "label": coefficient.name,
                    "value": coefficient.estimate,
                },
                {
                    "row_kind": "coefficient_standard_error",
                    "label": coefficient.name,
                    "value": coefficient.standard_error,
                },
                {
                    "row_kind": "coefficient_p_value",
                    "label": coefficient.name,
                    "value": coefficient.p_value,
                },
            ]
        )
    for matrix_row in model_matrix.rows:
        for column_name, value in matrix_row.encoded_values.items():
            rows.append(
                {
                    "row_kind": "model_matrix",
                    "label": f"{matrix_row.taxon}:{column_name}",
                    "value": value,
                }
            )
    return sorted(rows, key=lambda row: (str(row["row_kind"]), str(row["label"])))


def _phylogenetic_residual_parity_rows(report) -> list[dict[str, object]]:
    rows = [
        {
            "row_kind": "coefficient_estimate",
            "label": row.name,
            "value": row.estimate,
        }
        for row in report.coefficient_rows
    ]
    rows.extend(
        {
            "row_kind": "taxon_value",
            "label": row.taxon,
            "observed_value": row.observed_value,
            "fitted_value": row.fitted_value,
            "residual": row.residual,
        }
        for row in report.taxon_rows
    )
    return sorted(rows, key=lambda row: (str(row["row_kind"]), str(row["label"])))


def _phylogenetic_anova_parity_rows(report) -> list[dict[str, object]]:
    rows = [
        {
            "row_kind": "group_summary",
            "label": row.group,
            "taxon_count": row.taxon_count,
            "taxa": ",".join(row.taxa),
            "mean": row.mean,
            "variance": row.variance,
            "minimum": row.minimum,
            "maximum": row.maximum,
        }
        for row in report.group_rows
    ]
    rows.extend(
        {
            "row_kind": "pairwise_comparison",
            "label": f"{row.left_group}|{row.right_group}",
            "left_taxon_count": row.left_taxon_count,
            "right_taxon_count": row.right_taxon_count,
            "observed_t_statistic": row.observed_t_statistic,
            "uncorrected_p_value": row.uncorrected_p_value,
            "adjusted_p_value": row.adjusted_p_value,
        }
        for row in report.pairwise_rows
    )
    return sorted(rows, key=lambda row: (str(row["row_kind"]), str(row["label"])))


def build_comparative_case_payload(
    case: PhytoolsParityCase,
    *,
    tree_path: Path,
    traits_path: Path | None,
) -> tuple[dict[str, object], list[dict[str, object]] | None] | None:
    if case.operation == "comparative-pgls-brownian":
        if case.comparative_formula is None:
            raise ValueError("comparative-pgls-brownian requires comparative_formula")
        input_report = inspect_pgls_inputs(
            tree_path,
            traits_path,
            formula=case.comparative_formula,
            taxon_column=case.taxon_column,
        )
        model_matrix = build_pgls_model_matrix(
            tree_path,
            traits_path,
            formula=case.comparative_formula,
            taxon_column=case.taxon_column,
        )
        report = run_pgls(
            tree_path,
            traits_path,
            formula=case.comparative_formula,
            taxon_column=case.taxon_column,
            lambda_value=case.comparative_lambda_value or 1.0,
        )
        return (
            {
                "taxon_count": report.taxon_count,
                "trait_name": report.response,
                "formula": report.formula.formula,
                "analysis_taxon_count": len(report.taxa),
                "coefficient_count": len(report.coefficients),
                "model_matrix_row_count": len(model_matrix.rows),
                "model_matrix_column_count": len(model_matrix.encoded_columns),
                "categorical_predictor_count": len(input_report.categorical_predictors),
                "interaction_term_count": len(report.interaction_terms),
                "lambda_value": report.lambda_value,
                "lambda_estimation_mode": report.lambda_fit.mode,
                "log_likelihood": report.log_likelihood,
                "aic": report.aic,
                "residual_variance": report.residual_variance,
                "r_squared": report.r_squared,
                "diagnostic_outlier_count": len(report.diagnostics.outlier_taxa),
                "diagnostic_leverage_row_count": len(report.diagnostics.leverage_rows),
            },
            _pgls_parity_rows(model_matrix=model_matrix, report=report),
        )
    if case.operation == "phylogenetic-residuals":
        if case.comparative_predictors is None or len(case.comparative_predictors) != 1:
            raise ValueError(
                "phylogenetic-residuals requires one comparative predictor"
            )
        method = (
            "brownian"
            if case.comparative_lambda_value is not None
            and math.isclose(case.comparative_lambda_value, 1.0, abs_tol=1e-12)
            else "lambda"
        )
        report = summarize_phylogenetic_residuals(
            tree_path,
            traits_path,
            response=case.trait_name,
            predictor=case.comparative_predictors[0],
            taxon_column=case.taxon_column,
            method=method,
        )
        summary = {
            "taxon_count": report.analyzed_taxon_count,
            "trait_name": report.response,
            "predictor_name": report.predictor,
            "method": report.method,
            "excluded_taxon_count": len(report.excluded_taxa),
            "excluded_taxa": [row.taxon for row in report.excluded_taxa],
        }
        if method == "lambda":
            summary["lambda_value"] = report.lambda_value
            summary["log_likelihood"] = report.log_likelihood
        return summary, _phylogenetic_residual_parity_rows(report)
    if case.operation == "phylogenetic-anova":
        if case.comparative_predictors is None or len(case.comparative_predictors) != 1:
            raise ValueError("phylogenetic-anova requires one comparative predictor")
        report = summarize_phylogenetic_anova(
            tree_path,
            traits_path,
            response=case.trait_name,
            group=case.comparative_predictors[0],
            taxon_column=case.taxon_column,
            simulations=case.permutation_count or 199,
            seed=case.permutation_seed or 1,
        )
        return (
            {
                "taxon_count": report.analyzed_taxon_count,
                "trait_name": report.response,
                "group_column": report.group,
                "excluded_taxon_count": len(report.excluded_taxa),
                "excluded_taxa": [row.taxon for row in report.excluded_taxa],
                "group_count": report.group_count,
                "simulation_count": report.simulation_count,
                "seed": report.seed,
                "pairwise_adjustment_method": report.pairwise_adjustment_method,
                "brownian_sigma_squared": report.brownian_sigma_squared,
                "sum_of_squares_between": report.sum_of_squares_between,
                "sum_of_squares_within": report.sum_of_squares_within,
                "mean_square_between": report.mean_square_between,
                "mean_square_within": report.mean_square_within,
                "f_statistic": report.f_statistic,
                "p_value": report.p_value,
                "low_sample_group_count": report.low_sample_group_count,
            },
            _phylogenetic_anova_parity_rows(report),
        )
    return None
