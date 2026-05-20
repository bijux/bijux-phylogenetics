from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.pgls.brownian_covariance import (
    summarize_brownian_covariance_pgls,
    write_brownian_covariance_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_brownian_pgls_commands(comparative_subparsers: Any) -> None:
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


def run_brownian_pgls_command(args: Any) -> int | None:
    if args.comparative_command != "brownian-pgls":
        return None

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
