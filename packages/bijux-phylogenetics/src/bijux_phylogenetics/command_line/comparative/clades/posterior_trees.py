from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.pgls.posterior_tree import (
    run_posterior_tree_pgls,
    write_posterior_tree_pgls_coefficient_table,
    write_posterior_tree_pgls_summary_table,
    write_posterior_tree_pgls_tree_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_posterior_tree_regression_commands(comparative_subparsers: Any) -> None:
    comparative_posterior_pgls = comparative_subparsers.add_parser(
        "posterior-pgls",
        help="Fit one continuous-trait PGLS model across a posterior or bootstrap tree set.",
    )
    comparative_posterior_pgls.add_argument("tree_set", type=Path)
    comparative_posterior_pgls.add_argument("table", type=Path)
    comparative_posterior_pgls.add_argument("--response")
    comparative_posterior_pgls.add_argument("--predictors", nargs="+")
    comparative_posterior_pgls.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass + habitat'.",
    )
    comparative_posterior_pgls.add_argument("--taxon-column")
    comparative_posterior_pgls.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1 for each retained tree fit.",
    )
    comparative_posterior_pgls.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this leading fraction of the tree set before refitting.",
    )
    comparative_posterior_pgls.add_argument(
        "--significance-threshold",
        type=float,
        default=0.05,
        help="Treat coefficient p-values at or below this threshold as supported.",
    )
    comparative_posterior_pgls.add_argument(
        "--trees-out",
        type=Path,
        help="Write the per-tree posterior PGLS fit ledger as TSV or CSV.",
    )
    comparative_posterior_pgls.add_argument(
        "--coefficients-out",
        type=Path,
        help="Write the per-tree coefficient ledger as TSV or CSV.",
    )
    comparative_posterior_pgls.add_argument(
        "--summary-out",
        type=Path,
        help="Write the coefficient-distribution summary ledger as TSV or CSV.",
    )
    comparative_posterior_pgls.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior-tree PGLS report as JSON.",
    )
    _add_manifest_argument(comparative_posterior_pgls)


def run_posterior_tree_regression_command(
    args: Any,
    *,
    lambda_value: float | str,
) -> int | None:
    if args.comparative_command != "posterior-pgls":
        return None

    report = run_posterior_tree_pgls(
        args.tree_set,
        args.table,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
        burnin_fraction=args.burnin_fraction,
        significance_threshold=args.significance_threshold,
    )
    outputs: list[Path | str] = []
    if args.trees_out is not None:
        outputs.append(
            write_posterior_tree_pgls_tree_table(
                args.trees_out,
                report,
            )
        )
    if args.coefficients_out is not None:
        outputs.append(
            write_posterior_tree_pgls_coefficient_table(
                args.coefficients_out,
                report,
            )
        )
    if args.summary_out is not None:
        outputs.append(
            write_posterior_tree_pgls_summary_table(
                args.summary_out,
                report,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="comparative",
        inputs=[args.tree_set, args.table],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree_set, args.table],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "total_tree_count": report.total_tree_count,
                "burnin_tree_count": report.burnin_tree_count,
                "kept_tree_count": report.kept_tree_count,
                "analysis_taxon_count": len(report.analysis_taxa),
                "rooted_topology_count": report.rooted_topology_count,
                "unrooted_topology_count": report.unrooted_topology_count,
                "tree_fit_row_count": len(report.tree_rows),
                "coefficient_row_count": len(report.coefficient_rows),
                "coefficient_summary_count": len(report.coefficient_summaries),
                "stable_supported_term_count": sum(
                    row.conclusion_stability == "stable_supported"
                    for row in report.coefficient_summaries
                ),
                "direction_conflict_term_count": sum(
                    row.conclusion_stability == "direction_conflict"
                    for row in report.coefficient_summaries
                ),
                "lambda_mode": report.lambda_mode,
                "significance_threshold": report.significance_threshold,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
