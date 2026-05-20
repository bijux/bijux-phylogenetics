from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.pgls import run_pgls_multiple_testing
from bijux_phylogenetics.comparative.pgls.brownian_covariance import (
    summarize_brownian_covariance_pgls,
    write_brownian_covariance_table,
)
from bijux_phylogenetics.comparative.pgls.ou_covariance import (
    summarize_ou_covariance_pgls,
    write_ou_alpha_profile_table,
    write_ou_covariance_table,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .covariance_audit import (
    add_covariance_audit_pgls_commands,
    run_covariance_audit_pgls_command,
)
from .estimated_lambda import (
    add_estimated_lambda_pgls_commands,
    run_estimated_lambda_pgls_command,
)


def add_comparative_pgls_commands(comparative_subparsers: Any) -> None:
    add_covariance_audit_pgls_commands(comparative_subparsers)
    add_estimated_lambda_pgls_commands(comparative_subparsers)

    comparative_brownian_pgls = comparative_subparsers.add_parser(
        "brownian-pgls",
        help="Fit a PGLS model under fixed Brownian shared-path covariance.",
    )
    comparative_brownian_pgls.add_argument("tree", type=Path)
    comparative_brownian_pgls.add_argument("table", type=Path)
    comparative_brownian_pgls.add_argument("--response")
    comparative_brownian_pgls.add_argument("--predictors", nargs="+")
    comparative_brownian_pgls.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_brownian_pgls.add_argument("--taxon-column")
    comparative_brownian_pgls.add_argument(
        "--covariance-out",
        type=Path,
        help="Write the pairwise Brownian covariance ledger as TSV or CSV.",
    )
    comparative_brownian_pgls.add_argument(
        "--json", action="store_true", help="Emit the Brownian PGLS result as JSON."
    )
    _add_manifest_argument(comparative_brownian_pgls)

    comparative_ou_pgls = comparative_subparsers.add_parser(
        "ou-pgls",
        help="Fit a PGLS model under stationary-root OU covariance.",
    )
    comparative_ou_pgls.add_argument("tree", type=Path)
    comparative_ou_pgls.add_argument("table", type=Path)
    comparative_ou_pgls.add_argument("--response")
    comparative_ou_pgls.add_argument("--predictors", nargs="+")
    comparative_ou_pgls.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_ou_pgls.add_argument("--taxon-column")
    comparative_ou_pgls.add_argument(
        "--alpha",
        default="estimate",
        help="Use 'estimate' or a positive numeric OU alpha value.",
    )
    comparative_ou_pgls.add_argument(
        "--covariance-out",
        type=Path,
        help="Write the pairwise OU covariance ledger as TSV or CSV.",
    )
    comparative_ou_pgls.add_argument(
        "--alpha-profile-out",
        type=Path,
        help="Write the fitted OU alpha likelihood profile as TSV or CSV.",
    )
    comparative_ou_pgls.add_argument(
        "--json",
        action="store_true",
        help="Emit the OU covariance PGLS result as JSON.",
    )
    _add_manifest_argument(comparative_ou_pgls)

    comparative_multiple_testing = comparative_subparsers.add_parser(
        "multiple-testing",
        help="Adjust PGLS coefficient p-values across many response traits.",
    )
    comparative_multiple_testing.add_argument("tree", type=Path)
    comparative_multiple_testing.add_argument("table", type=Path)
    comparative_multiple_testing.add_argument("--responses", nargs="+", required=True)
    comparative_multiple_testing.add_argument("--predictors", nargs="+", required=True)
    comparative_multiple_testing.add_argument("--taxon-column")
    comparative_multiple_testing.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_multiple_testing.add_argument(
        "--json", action="store_true", help="Emit the correction report as JSON."
    )
    _add_manifest_argument(comparative_multiple_testing)


def run_comparative_pgls_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    del parser
    lambda_value: float | str
    if hasattr(args, "lambda_value"):
        if args.lambda_value == "estimate":
            lambda_value = "estimate"
        else:
            lambda_value = float(args.lambda_value)
    else:
        lambda_value = "estimate"
    covariance_audit_result = run_covariance_audit_pgls_command(
        args,
        lambda_value=lambda_value,
    )
    if covariance_audit_result is not None:
        return covariance_audit_result

    estimated_lambda_result = run_estimated_lambda_pgls_command(
        args,
        lambda_value=lambda_value,
    )
    if estimated_lambda_result is not None:
        return estimated_lambda_result
    if args.comparative_command == "multiple-testing":
        report = run_pgls_multiple_testing(
            args.tree,
            args.table,
            responses=list(args.responses),
            predictors=list(args.predictors),
            taxon_column=args.taxon_column,
            lambda_value=lambda_value,
        )
        outputs = _finalize_outputs(
            args,
            command="comparative",
            inputs=[args.tree, args.table],
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "response_count": len(report.responses),
                    "test_count": len(report.rows),
                    "family_size": report.family_size,
                    "raw_significant_count": report.raw_significant_count,
                    "significant_count": sum(1 for row in report.rows if row.significant),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "brownian-pgls":
        report = summarize_brownian_covariance_pgls(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
        )
        outputs: list[Path | str] = []
        if args.covariance_out is not None:
            outputs.append(
                write_brownian_covariance_table(
                    args.covariance_out,
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
                    "taxon_count": report.taxon_count,
                    "predictor_count": len(report.model.predictors),
                    "coefficient_count": len(report.model.coefficients),
                    "covariance_row_count": len(report.rows),
                    "lambda_value": report.model.lambda_value,
                    "covariance_model": "brownian-shared-path",
                    "tree_is_ultrametric": report.tree_is_ultrametric,
                    "minimum_root_to_tip_depth": report.minimum_root_to_tip_depth,
                    "maximum_root_to_tip_depth": report.maximum_root_to_tip_depth,
                    "raw_log_determinant": report.raw_log_determinant,
                    "positive_definite_before_stabilization": (
                        report.positive_definite_before_stabilization
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "ou-pgls":
        report = summarize_ou_covariance_pgls(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
            alpha=args.alpha,
        )
        outputs: list[Path | str] = []
        if args.covariance_out is not None:
            outputs.append(write_ou_covariance_table(args.covariance_out, report))
        if args.alpha_profile_out is not None:
            outputs.append(
                write_ou_alpha_profile_table(
                    args.alpha_profile_out,
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
                    "taxon_count": report.taxon_count,
                    "predictor_count": len(report.model.predictors),
                    "coefficient_count": len(report.model.coefficients),
                    "covariance_row_count": len(report.rows),
                    "alpha": report.alpha,
                    "alpha_estimation_mode": report.alpha_estimation_mode,
                    "alpha_profile_point_count": len(report.alpha_profile_rows),
                    "alpha_lower_95_confidence_interval": (
                        report.lower_95_confidence_interval
                    ),
                    "alpha_upper_95_confidence_interval": (
                        report.upper_95_confidence_interval
                    ),
                    "covariance_model": "ou-stationary-root",
                    "tree_is_ultrametric": report.tree_is_ultrametric,
                    "raw_log_determinant": report.raw_log_determinant,
                    "positive_definite_before_stabilization": (
                        report.positive_definite_before_stabilization
                    ),
                    "log_likelihood": report.model.log_likelihood,
                    "aic": report.model.aic,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    return None
