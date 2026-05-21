from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.pgls import (
    inspect_pgls_inputs,
    run_pgls,
    write_pgls_model_matrix_table,
)
from bijux_phylogenetics.comparative.pgls.categorical_contrasts import (
    summarize_pgls_categorical_contrasts,
    write_pgls_categorical_contrast_table,
)
from bijux_phylogenetics.comparative.pgls.interaction_coefficients import (
    summarize_pgls_interaction_coefficients,
    write_pgls_interaction_coefficient_table,
)
from bijux_phylogenetics.comparative.pgls.lambda_fit import (
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_estimated_lambda_pgls_commands(comparative_subparsers: Any) -> None:
    comparative_pgls = comparative_subparsers.add_parser(
        "pgls",
        help="Fit a phylogenetic generalized least-squares model.",
    )
    comparative_pgls.add_argument("tree", type=Path)
    comparative_pgls.add_argument("table", type=Path)
    comparative_pgls.add_argument("--response")
    comparative_pgls.add_argument("--predictors", nargs="+")
    comparative_pgls.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_pgls.add_argument("--taxon-column")
    comparative_pgls.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_pgls.add_argument(
        "--model-matrix-out",
        type=Path,
        help="Write the encoded comparative model matrix as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--categorical-contrasts-out",
        type=Path,
        help="Write categorical predictor contrast rows as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--interaction-coefficients-out",
        type=Path,
        help="Write interaction coefficient rows as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--lambda-profile-out",
        type=Path,
        help="Write the fitted Pagel lambda likelihood profile as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--json", action="store_true", help="Emit the model result as JSON."
    )
    _add_manifest_argument(comparative_pgls)


def run_estimated_lambda_pgls_command(
    args: Any,
    *,
    lambda_value: float | str,
) -> int | None:
    if args.comparative_command != "pgls":
        return None

    input_report = inspect_pgls_inputs(
        args.tree,
        args.table,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        taxon_column=args.taxon_column,
    )
    report = run_pgls(
        args.tree,
        args.table,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
    )
    categorical_contrasts = summarize_pgls_categorical_contrasts(
        args.tree,
        args.table,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
    )
    interaction_coefficients = summarize_pgls_interaction_coefficients(
        args.tree,
        args.table,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
    )
    outputs: list[Path | str] = []
    if args.model_matrix_out is not None:
        outputs.append(
            write_pgls_model_matrix_table(
                args.model_matrix_out,
                input_report.model_matrix,
            )
        )
    if args.categorical_contrasts_out is not None:
        outputs.append(
            write_pgls_categorical_contrast_table(
                args.categorical_contrasts_out,
                categorical_contrasts,
            )
        )
    if args.interaction_coefficients_out is not None:
        outputs.append(
            write_pgls_interaction_coefficient_table(
                args.interaction_coefficients_out,
                interaction_coefficients,
            )
        )
    if args.lambda_profile_out is not None:
        outputs.append(
            write_pgls_lambda_profile_table(
                args.lambda_profile_out,
                report.lambda_fit,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="comparative",
        inputs=[args.tree, args.table],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
            warnings=input_report.warnings,
            metrics={
                "taxon_count": report.taxon_count,
                "predictor_count": len(report.predictors),
                "coefficient_count": len(report.coefficients),
                "confidence_interval_count": len(report.coefficients),
                "categorical_contrast_predictor_count": (
                    categorical_contrasts.categorical_predictor_count
                ),
                "categorical_contrast_row_count": len(categorical_contrasts.rows),
                "interaction_term_count": (
                    interaction_coefficients.interaction_term_count
                ),
                "interaction_coefficient_row_count": len(interaction_coefficients.rows),
                "intercept_included": input_report.formula.include_intercept,
                "model_matrix_row_count": input_report.model_matrix.row_count,
                "model_matrix_column_count": len(
                    input_report.model_matrix.encoded_columns
                ),
                "residual_degrees_of_freedom": (
                    report.coefficients[0].degrees_of_freedom
                    if report.coefficients
                    else 0
                ),
                "coefficient_inference_distribution": (
                    report.coefficients[0].inference_distribution
                    if report.coefficients
                    else None
                ),
                "encoded_predictor_count": len(
                    input_report.model_matrix.encoded_columns
                )
                - (1 if input_report.formula.include_intercept else 0),
                "categorical_predictor_count": len(input_report.categorical_predictors),
                "transformed_term_count": len(
                    input_report.formula_audit.transformed_terms
                ),
                "lambda_value": report.lambda_value,
                "lambda_estimation_mode": report.lambda_fit.mode,
                "lambda_profile_point_count": len(report.lambda_fit.profile_rows),
                "lambda_lower_95_confidence_interval": (
                    report.lambda_fit.lower_95_confidence_interval
                ),
                "lambda_upper_95_confidence_interval": (
                    report.lambda_fit.upper_95_confidence_interval
                ),
                "aic": report.aic,
                "r_squared": report.r_squared,
            },
            data={
                "inputs": input_report,
                "model": report,
                "categorical_contrasts": categorical_contrasts,
                "interaction_coefficients": interaction_coefficients,
            },
        ),
        json_output=args.json,
    )
    return 0
