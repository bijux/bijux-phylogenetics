from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.clade_residuals import (
    analyze_comparative_residual_clades,
    write_comparative_residual_clade_table,
    write_comparative_residual_taxon_table,
)
from bijux_phylogenetics.comparative.clade_stability import (
    analyze_comparative_clade_stability,
    write_comparative_clade_coefficient_change_table,
    write_comparative_clade_stability_table,
)
from bijux_phylogenetics.comparative.pgls.posterior_tree import (
    run_posterior_tree_pgls,
    write_posterior_tree_pgls_coefficient_table,
    write_posterior_tree_pgls_summary_table,
    write_posterior_tree_pgls_tree_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_comparative_support_commands(comparative_subparsers: Any) -> None:
    comparative_clade_residuals = comparative_subparsers.add_parser(
        "clade-residuals",
        help="Aggregate comparative model residuals across internal clades.",
    )
    comparative_clade_residuals.add_argument("tree", type=Path)
    comparative_clade_residuals.add_argument("table", type=Path)
    comparative_clade_residuals.add_argument("--response")
    comparative_clade_residuals.add_argument("--predictors", nargs="+")
    comparative_clade_residuals.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass + habitat'.",
    )
    comparative_clade_residuals.add_argument("--taxon-column")
    comparative_clade_residuals.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1. Binary-response residual aggregation requires a numeric value.",
    )
    comparative_clade_residuals.add_argument(
        "--taxa-out",
        type=Path,
        help="Write the analyzed taxon residual ledger as TSV or CSV.",
    )
    comparative_clade_residuals.add_argument(
        "--clades-out",
        type=Path,
        help="Write the internal clade residual aggregation ledger as TSV or CSV.",
    )
    comparative_clade_residuals.add_argument(
        "--json", action="store_true", help="Emit the clade residual report as JSON."
    )
    _add_manifest_argument(comparative_clade_residuals)

    comparative_clade_stability = comparative_subparsers.add_parser(
        "clade-stability",
        help="Refit one comparative model after removing each major internal clade.",
    )
    comparative_clade_stability.add_argument("tree", type=Path)
    comparative_clade_stability.add_argument("table", type=Path)
    comparative_clade_stability.add_argument("--response")
    comparative_clade_stability.add_argument("--predictors", nargs="+")
    comparative_clade_stability.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass + habitat'.",
    )
    comparative_clade_stability.add_argument("--taxon-column")
    comparative_clade_stability.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1. Binary-response clade stability requires a numeric value.",
    )
    comparative_clade_stability.add_argument(
        "--clades-out",
        type=Path,
        help="Write the leave-one-clade-out stability summary ledger as TSV or CSV.",
    )
    comparative_clade_stability.add_argument(
        "--terms-out",
        type=Path,
        help="Write the coefficient-delta ledger across clade removals as TSV or CSV.",
    )
    comparative_clade_stability.add_argument(
        "--json", action="store_true", help="Emit the clade-stability report as JSON."
    )
    _add_manifest_argument(comparative_clade_stability)

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


def run_comparative_support_command(
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
    if args.comparative_command == "clade-residuals":
        report = analyze_comparative_residual_clades(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
            lambda_value=lambda_value,
        )
        outputs: list[Path | str] = []
        if args.taxa_out is not None:
            outputs.append(
                write_comparative_residual_taxon_table(
                    args.taxa_out,
                    report,
                )
            )
        if args.clades_out is not None:
            outputs.append(
                write_comparative_residual_clade_table(
                    args.clades_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
        )
        top_clade = (
            min(report.clade_rows, key=lambda row: row.rank).clade_id
            if report.clade_rows
            else None
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "model_family": report.model_family,
                    "taxon_count": len(report.taxon_rows),
                    "clade_count": len(report.clade_rows),
                    "residual_heavy_clade_count": len(report.residual_heavy_clades),
                    "top_influential_clade": top_clade,
                    "standardized_residual_method": report.standardized_residual_method,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "clade-stability":
        report = analyze_comparative_clade_stability(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
            lambda_value=lambda_value,
        )
        outputs: list[Path | str] = []
        if args.clades_out is not None:
            outputs.append(
                write_comparative_clade_stability_table(
                    args.clades_out,
                    report,
                )
            )
        if args.terms_out is not None:
            outputs.append(
                write_comparative_clade_coefficient_change_table(
                    args.terms_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
        )
        top_clade = (
            min(
                (
                    row
                    for row in report.clade_rows
                    if row.fit_status == "fit" and row.rank > 0
                ),
                key=lambda row: row.rank,
            ).clade_id
            if any(
                row.fit_status == "fit" and row.rank > 0
                for row in report.clade_rows
            )
            else None
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "model_family": report.model_family,
                    "baseline_taxon_count": len(report.baseline_taxa),
                    "baseline_term_count": report.baseline_term_count,
                    "candidate_clade_count": report.candidate_clade_count,
                    "blocked_clade_count": report.blocked_clade_count,
                    "coefficient_change_row_count": len(report.coefficient_rows),
                    "top_influential_clade": top_clade,
                    "major_clade_fraction": report.major_clade_fraction,
                    "minimum_major_clade_size": report.minimum_major_clade_size,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "posterior-pgls":
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
    return None
