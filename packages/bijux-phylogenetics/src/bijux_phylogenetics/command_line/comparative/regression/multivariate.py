from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.regression import (
    run_multivariate_comparative_regression,
    write_multivariate_excluded_taxa_table,
    write_multivariate_residual_association_table,
    write_multivariate_residual_correlation_table,
    write_multivariate_residual_covariance_table,
    write_multivariate_response_coefficient_table,
    write_multivariate_response_model_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_multivariate_regression_command(comparative_subparsers: Any) -> None:
    comparative_multivariate = comparative_subparsers.add_parser(
        "multivariate",
        help="Fit shared-taxon comparative regressions across multiple response traits.",
    )
    comparative_multivariate.add_argument("tree", type=Path)
    comparative_multivariate.add_argument("table", type=Path)
    comparative_multivariate.add_argument("--responses", nargs="+", required=True)
    comparative_multivariate.add_argument("--predictors", nargs="+", required=True)
    comparative_multivariate.add_argument("--taxon-column")
    comparative_multivariate.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_multivariate.add_argument(
        "--covariance-out",
        type=Path,
        help="Write the residual covariance ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--correlation-out",
        type=Path,
        help="Write the residual correlation ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--associations-out",
        type=Path,
        help="Write the residual trait-association ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--coefficients-out",
        type=Path,
        help="Write the per-response coefficient ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--response-models-out",
        type=Path,
        help="Write the per-response model summary ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write the explicit excluded-taxa ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--json",
        action="store_true",
        help="Emit the multivariate regression report as JSON.",
    )
    _add_manifest_argument(comparative_multivariate)


def run_multivariate_regression_command(
    args: Any,
    *,
    lambda_value: float | str,
) -> int | None:
    if args.comparative_command != "multivariate":
        return None

    report = run_multivariate_comparative_regression(
        args.tree,
        args.table,
        responses=list(args.responses),
        predictors=list(args.predictors),
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
    )
    outputs: list[Path | str] = []
    if args.covariance_out is not None:
        outputs.append(
            write_multivariate_residual_covariance_table(args.covariance_out, report)
        )
    if args.correlation_out is not None:
        outputs.append(
            write_multivariate_residual_correlation_table(
                args.correlation_out,
                report,
            )
        )
    if args.associations_out is not None:
        outputs.append(
            write_multivariate_residual_association_table(
                args.associations_out,
                report,
            )
        )
    if args.coefficients_out is not None:
        outputs.append(
            write_multivariate_response_coefficient_table(
                args.coefficients_out,
                report,
            )
        )
    if args.response_models_out is not None:
        outputs.append(
            write_multivariate_response_model_table(
                args.response_models_out,
                report,
            )
        )
    if args.excluded_taxa_out is not None:
        outputs.append(
            write_multivariate_excluded_taxa_table(
                args.excluded_taxa_out,
                report,
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
            metrics={
                "response_count": len(report.responses),
                "predictor_count": len(report.predictors),
                "analysis_taxa": len(report.analysis_taxa),
                "excluded_taxa": len(report.excluded_taxa),
                "residual_covariance_response_count": (
                    report.covariance_diagnostics.response_count
                ),
                "residual_covariance_matrix_rank": (
                    report.covariance_diagnostics.matrix_rank
                ),
                "residual_covariance_condition_number": (
                    report.covariance_diagnostics.condition_number
                ),
                "residual_covariance_singular": (
                    report.covariance_diagnostics.is_singular
                ),
                "residual_covariance_near_singular": (
                    report.covariance_diagnostics.is_near_singular
                ),
                "residual_covariance_row_count": len(report.covariance_rows),
                "residual_correlation_row_count": len(report.correlation_rows),
                "residual_association_count": len(report.association_rows),
                "response_model_count": len(report.response_model_rows),
                "coefficient_row_count": len(report.coefficient_rows),
                "warning_count": len(report.warnings),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
