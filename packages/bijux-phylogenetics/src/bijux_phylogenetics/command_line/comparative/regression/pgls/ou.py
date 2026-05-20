from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.pgls.ou_covariance import (
    summarize_ou_covariance_pgls,
    write_ou_alpha_profile_table,
    write_ou_covariance_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_ou_pgls_commands(comparative_subparsers: Any) -> None:
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


def run_ou_pgls_command(args: Any) -> int | None:
    if args.comparative_command != "ou-pgls":
        return None

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
