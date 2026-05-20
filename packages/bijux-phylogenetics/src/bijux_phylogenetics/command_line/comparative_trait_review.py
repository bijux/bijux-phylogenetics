from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.clades.traits import (
    summarize_clade_traits,
    write_clade_trait_clade_table,
    write_clade_trait_exclusion_table,
    write_clade_trait_summary_table,
)
from bijux_phylogenetics.comparative.phylogenetic_anova import (
    summarize_phylogenetic_anova,
    write_phylogenetic_anova_exclusion_table,
    write_phylogenetic_anova_group_table,
    write_phylogenetic_anova_pairwise_table,
    write_phylogenetic_anova_simulation_table,
    write_phylogenetic_anova_summary_table,
)
from bijux_phylogenetics.comparative.phylogenetic_residuals import (
    summarize_phylogenetic_residuals,
    write_phylogenetic_residual_coefficient_table,
    write_phylogenetic_residual_exclusion_table,
    write_phylogenetic_residual_summary_table,
    write_phylogenetic_residual_taxon_table,
)
from bijux_phylogenetics.comparative.traits.imputation import (
    summarize_trait_imputation,
    write_trait_imputation_exclusion_table,
    write_trait_imputation_holdout_table,
    write_trait_imputation_summary_table,
    write_trait_imputation_table,
)
from bijux_phylogenetics.comparative.traits.outliers import (
    summarize_trait_outliers,
    write_trait_outlier_exclusion_table,
    write_trait_outlier_summary_table,
    write_trait_outlier_taxon_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_comparative_trait_review_commands(comparative_subparsers: Any) -> None:
    comparative_clade_traits = comparative_subparsers.add_parser(
        "clade-traits",
        help="Summarize one continuous or categorical trait across internal clades.",
    )
    comparative_clade_traits.add_argument("tree", type=Path)
    comparative_clade_traits.add_argument("table", type=Path)
    comparative_clade_traits.add_argument("--trait", required=True)
    comparative_clade_traits.add_argument("--taxon-column")
    comparative_clade_traits.add_argument(
        "--trait-kind",
        choices=("auto", "continuous", "categorical"),
        default="auto",
        help="Infer trait kind automatically or force continuous/categorical handling.",
    )
    comparative_clade_traits.add_argument(
        "--min-clade-size",
        type=int,
        default=2,
        help="Only summarize internal clades with at least this many analyzed taxa.",
    )
    comparative_clade_traits.add_argument(
        "--summary-out",
        type=Path,
        help="Write one clade-trait summary ledger as TSV or CSV.",
    )
    comparative_clade_traits.add_argument(
        "--clades-out",
        type=Path,
        help="Write one internal clade-trait ledger as TSV or CSV.",
    )
    comparative_clade_traits.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for clade trait summarization as TSV or CSV.",
    )
    comparative_clade_traits.add_argument(
        "--json", action="store_true", help="Emit the clade-trait report as JSON."
    )
    _add_manifest_argument(comparative_clade_traits)
    comparative_phylogenetic_residuals = comparative_subparsers.add_parser(
        "phylogenetic-residuals",
        help="Summarize tree-aware fitted values and residuals for one continuous response and predictor.",
    )
    comparative_phylogenetic_residuals.add_argument("tree", type=Path)
    comparative_phylogenetic_residuals.add_argument("table", type=Path)
    comparative_phylogenetic_residuals.add_argument("--response", required=True)
    comparative_phylogenetic_residuals.add_argument("--predictor", required=True)
    comparative_phylogenetic_residuals.add_argument("--taxon-column")
    comparative_phylogenetic_residuals.add_argument(
        "--method",
        choices=("brownian", "lambda"),
        default="lambda",
        help="Use fixed Brownian covariance or estimate Pagel lambda before computing residuals.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--summary-out",
        type=Path,
        help="Write one phylogenetic-residual summary ledger as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--residuals-out",
        type=Path,
        help="Write one taxon-level phylogenetic-residual ledger as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--coefficients-out",
        type=Path,
        help="Write one phylogenetic-residual coefficient ledger as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for phylogenetic residual review as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--json",
        action="store_true",
        help="Emit the phylogenetic-residual review as JSON.",
    )
    _add_manifest_argument(comparative_phylogenetic_residuals)
    comparative_phylogenetic_anova = comparative_subparsers.add_parser(
        "phylogenetic-anova",
        help="Run a simulation-based phylogenetic ANOVA over one continuous response and one categorical group.",
    )
    comparative_phylogenetic_anova.add_argument("tree", type=Path)
    comparative_phylogenetic_anova.add_argument("table", type=Path)
    comparative_phylogenetic_anova.add_argument("--response", required=True)
    comparative_phylogenetic_anova.add_argument("--group", required=True)
    comparative_phylogenetic_anova.add_argument("--taxon-column")
    comparative_phylogenetic_anova.add_argument(
        "--simulations",
        type=int,
        default=1000,
        help="Number of observed-plus-null F statistics to evaluate.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Seed for the Brownian null simulation sequence.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--summary-out",
        type=Path,
        help="Write one phylogenetic-ANOVA summary ledger as TSV or CSV.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--groups-out",
        type=Path,
        help="Write one group-summary ledger for phylogenetic ANOVA as TSV or CSV.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--pairwise-out",
        type=Path,
        help="Write one pairwise-comparison ledger for phylogenetic ANOVA as TSV or CSV.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--simulations-out",
        type=Path,
        help="Write one observed-plus-null F-statistic ledger for phylogenetic ANOVA as TSV or CSV.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for phylogenetic ANOVA as TSV or CSV.",
    )
    comparative_phylogenetic_anova.add_argument(
        "--json",
        action="store_true",
        help="Emit the phylogenetic-ANOVA report as JSON.",
    )
    _add_manifest_argument(comparative_phylogenetic_anova)
    comparative_trait_outliers = comparative_subparsers.add_parser(
        "trait-outliers",
        help="Rank continuous-trait taxa by conditional phylogenetic residual size.",
    )
    comparative_trait_outliers.add_argument("tree", type=Path)
    comparative_trait_outliers.add_argument("table", type=Path)
    comparative_trait_outliers.add_argument("--trait", required=True)
    comparative_trait_outliers.add_argument("--taxon-column")
    comparative_trait_outliers.add_argument(
        "--summary-out",
        type=Path,
        help="Write one trait-outlier summary ledger as TSV or CSV.",
    )
    comparative_trait_outliers.add_argument(
        "--outliers-out",
        type=Path,
        help="Write one ranked taxon outlier ledger as TSV or CSV.",
    )
    comparative_trait_outliers.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for trait outlier review as TSV or CSV.",
    )
    comparative_trait_outliers.add_argument(
        "--json",
        action="store_true",
        help="Emit the trait outlier review as JSON.",
    )
    _add_manifest_argument(comparative_trait_outliers)
    comparative_trait_imputation = comparative_subparsers.add_parser(
        "trait-imputation",
        help="Impute missing continuous-trait values under a Brownian phylogenetic model.",
    )
    comparative_trait_imputation.add_argument("tree", type=Path)
    comparative_trait_imputation.add_argument("table", type=Path)
    comparative_trait_imputation.add_argument("--trait", required=True)
    comparative_trait_imputation.add_argument("--taxon-column")
    comparative_trait_imputation.add_argument(
        "--summary-out",
        type=Path,
        help="Write one Brownian trait-imputation summary ledger as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--imputations-out",
        type=Path,
        help="Write one imputed-value ledger as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--holdout-out",
        type=Path,
        help="Write one leave-one-observed-out validation ledger as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for trait imputation as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--json",
        action="store_true",
        help="Emit the trait-imputation review as JSON.",
    )
    _add_manifest_argument(comparative_trait_imputation)


def run_comparative_trait_review_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    if args.comparative_command == "clade-traits":
        report = summarize_clade_traits(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            minimum_clade_size=args.min_clade_size,
            trait_kind=args.trait_kind,
        )
        if args.summary_out:
            write_clade_trait_summary_table(args.summary_out, report)
        if args.clades_out:
            write_clade_trait_clade_table(args.clades_out, report)
        if args.excluded_taxa_out:
            write_clade_trait_exclusion_table(
                args.excluded_taxa_out,
                report,
            )
        outputs = _finalize_outputs(
            args, command="comparative", inputs=[args.tree, args.table]
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "trait_kind": report.trait_kind,
                    "clade_count": len(report.clade_rows),
                    "exceptional_clade_count": len(report.exceptional_clades),
                    "top_exceptional_clade": report.top_exceptional_clade,
                    "top_exceptionality_score": report.top_exceptionality_score,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "phylogenetic-residuals":
        report = summarize_phylogenetic_residuals(
            args.tree,
            args.table,
            response=args.response,
            predictor=args.predictor,
            taxon_column=args.taxon_column,
            method=args.method,
        )
        if args.summary_out:
            write_phylogenetic_residual_summary_table(args.summary_out, report)
        if args.residuals_out:
            write_phylogenetic_residual_taxon_table(args.residuals_out, report)
        if args.coefficients_out:
            write_phylogenetic_residual_coefficient_table(
                args.coefficients_out,
                report,
            )
        if args.excluded_taxa_out:
            write_phylogenetic_residual_exclusion_table(
                args.excluded_taxa_out,
                report,
            )
        outputs = _finalize_outputs(
            args, command="comparative", inputs=[args.tree, args.table]
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "method": report.method,
                    "outlier_count": len(report.outlier_taxa),
                    "top_outlier_taxon": (
                        None
                        if not report.taxon_rows
                        else max(
                            report.taxon_rows,
                            key=lambda row: row.abs_standardized_residual,
                        ).taxon
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "phylogenetic-anova":
        report = summarize_phylogenetic_anova(
            args.tree,
            args.table,
            response=args.response,
            group=args.group,
            taxon_column=args.taxon_column,
            simulations=args.simulations,
            seed=args.seed,
        )
        if args.summary_out:
            write_phylogenetic_anova_summary_table(args.summary_out, report)
        if args.groups_out:
            write_phylogenetic_anova_group_table(args.groups_out, report)
        if args.pairwise_out:
            write_phylogenetic_anova_pairwise_table(args.pairwise_out, report)
        if args.simulations_out:
            write_phylogenetic_anova_simulation_table(
                args.simulations_out,
                report,
            )
        if args.excluded_taxa_out:
            write_phylogenetic_anova_exclusion_table(
                args.excluded_taxa_out,
                report,
            )
        outputs = _finalize_outputs(
            args, command="comparative", inputs=[args.tree, args.table]
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "group_count": report.group_count,
                    "simulation_count": report.simulation_count,
                    "p_value": report.p_value,
                    "low_sample_group_count": report.low_sample_group_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "trait-outliers":
        report = summarize_trait_outliers(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
        )
        if args.summary_out:
            write_trait_outlier_summary_table(args.summary_out, report)
        if args.outliers_out:
            write_trait_outlier_taxon_table(args.outliers_out, report)
        if args.excluded_taxa_out:
            write_trait_outlier_exclusion_table(
                args.excluded_taxa_out,
                report,
            )
        outputs = _finalize_outputs(
            args, command="comparative", inputs=[args.tree, args.table]
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "selected_model": report.selected_model,
                    "outlier_count": len(report.outlier_taxa),
                    "top_outlier_taxon": report.top_outlier_taxon,
                    "top_abs_standardized_residual": (
                        report.top_abs_standardized_residual
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "trait-imputation":
        report = summarize_trait_imputation(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
        )
        if args.summary_out:
            write_trait_imputation_summary_table(args.summary_out, report)
        if args.imputations_out:
            write_trait_imputation_table(args.imputations_out, report)
        if args.holdout_out:
            write_trait_imputation_holdout_table(args.holdout_out, report)
        if args.excluded_taxa_out:
            write_trait_imputation_exclusion_table(
                args.excluded_taxa_out,
                report,
            )
        outputs = _finalize_outputs(
            args, command="comparative", inputs=[args.tree, args.table]
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "observed_taxon_count": report.observed_taxon_count,
                    "imputed_taxon_count": len(report.imputation_rows),
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "holdout_validation_status": report.holdout_validation_status,
                    "holdout_count": len(report.holdout_rows),
                    "holdout_mean_absolute_error": (
                        report.holdout_mean_absolute_error
                    ),
                    "holdout_interval_coverage": (
                        report.holdout_interval_coverage
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
