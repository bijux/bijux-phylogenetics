from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.continuous import (
    summarize_brownian_trait_evolution,
    summarize_early_burst_trait_evolution,
    summarize_ou_trait_evolution,
    write_brownian_trait_evolution_exclusion_table,
    write_brownian_trait_evolution_summary_table,
    write_early_burst_rate_change_profile_table,
    write_early_burst_trait_evolution_comparison_table,
    write_early_burst_trait_evolution_exclusion_table,
    write_early_burst_trait_evolution_summary_table,
    write_ou_trait_evolution_exclusion_table,
    write_ou_trait_evolution_summary_table,
)
from bijux_phylogenetics.comparative.traits.rate_through_time import (
    summarize_trait_rate_through_time,
    write_trait_rate_through_time_exclusion_table,
    write_trait_rate_through_time_interval_table,
    write_trait_rate_through_time_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_continuous_model_comparative_evolution_commands(
    comparative_subparsers: Any,
) -> None:
    comparative_brownian = comparative_subparsers.add_parser(
        "brownian",
        help="Fit a standalone Brownian-motion continuous-trait model.",
    )
    comparative_brownian.add_argument("tree", type=Path)
    comparative_brownian.add_argument("table", type=Path)
    comparative_brownian.add_argument("--trait", required=True)
    comparative_brownian.add_argument("--taxon-column")
    comparative_brownian.add_argument(
        "--summary-out",
        type=Path,
        help="Write one Brownian trait-evolution summary ledger as TSV or CSV.",
    )
    comparative_brownian.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for the Brownian trait fit as TSV or CSV.",
    )
    comparative_brownian.add_argument(
        "--json", action="store_true", help="Emit the Brownian model fit as JSON."
    )
    _add_manifest_argument(comparative_brownian)

    comparative_ou = comparative_subparsers.add_parser(
        "ou",
        help="Fit a standalone Ornstein-Uhlenbeck continuous-trait model.",
    )
    comparative_ou.add_argument("tree", type=Path)
    comparative_ou.add_argument("table", type=Path)
    comparative_ou.add_argument("--trait", required=True)
    comparative_ou.add_argument("--taxon-column")
    comparative_ou.add_argument(
        "--summary-out",
        type=Path,
        help="Write one OU trait-evolution summary ledger as TSV or CSV.",
    )
    comparative_ou.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for the OU trait fit as TSV or CSV.",
    )
    comparative_ou.add_argument(
        "--json", action="store_true", help="Emit the OU model fit as JSON."
    )
    _add_manifest_argument(comparative_ou)

    comparative_early_burst = comparative_subparsers.add_parser(
        "early-burst",
        help="Fit a standalone early-burst continuous-trait model with BM/OU comparison.",
    )
    comparative_early_burst.add_argument("tree", type=Path)
    comparative_early_burst.add_argument("table", type=Path)
    comparative_early_burst.add_argument("--trait", required=True)
    comparative_early_burst.add_argument("--taxon-column")
    comparative_early_burst.add_argument(
        "--summary-out",
        type=Path,
        help="Write one early-burst trait-evolution summary ledger as TSV or CSV.",
    )
    comparative_early_burst.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for the early-burst trait fit as TSV or CSV.",
    )
    comparative_early_burst.add_argument(
        "--comparison-out",
        type=Path,
        help="Write one BM/OU/early-burst comparison ledger as TSV or CSV.",
    )
    comparative_early_burst.add_argument(
        "--profile-out",
        type=Path,
        help="Write one bounded rate-change likelihood profile as TSV or CSV.",
    )
    comparative_early_burst.add_argument(
        "--json", action="store_true", help="Emit the early-burst model fit as JSON."
    )
    _add_manifest_argument(comparative_early_burst)

    comparative_rate_through_time = comparative_subparsers.add_parser(
        "rate-through-time",
        help="Summarize how trait-rate evidence changes across tree depth intervals.",
    )
    comparative_rate_through_time.add_argument("tree", type=Path)
    comparative_rate_through_time.add_argument("table", type=Path)
    comparative_rate_through_time.add_argument("--trait", required=True)
    comparative_rate_through_time.add_argument("--taxon-column")
    comparative_rate_through_time.add_argument(
        "--interval-count",
        type=int,
        default=5,
        help="Number of equal-width depth intervals used for the rate-through-time ledger.",
    )
    comparative_rate_through_time.add_argument(
        "--summary-out",
        type=Path,
        help="Write one rate-through-time summary ledger as TSV or CSV.",
    )
    comparative_rate_through_time.add_argument(
        "--intervals-out",
        type=Path,
        help="Write one interval rate ledger as TSV or CSV.",
    )
    comparative_rate_through_time.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for the rate-through-time fit as TSV or CSV.",
    )
    comparative_rate_through_time.add_argument(
        "--json",
        action="store_true",
        help="Emit the rate-through-time report as JSON.",
    )
    _add_manifest_argument(comparative_rate_through_time)


def run_continuous_model_comparative_evolution_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    del parser
    if args.comparative_command == "brownian":
        report = summarize_brownian_trait_evolution(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
        )
        if args.summary_out:
            write_brownian_trait_evolution_summary_table(args.summary_out, report)
        if args.excluded_taxa_out:
            write_brownian_trait_evolution_exclusion_table(
                args.excluded_taxa_out,
                report,
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
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "root_state": report.root_state,
                    "sigma_squared": report.sigma_squared,
                    "rate": report.sigma_squared,
                    "log_likelihood": report.log_likelihood,
                    "aic": report.aic,
                    "aicc": report.aicc,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command == "ou":
        report = summarize_ou_trait_evolution(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
        )
        if args.summary_out:
            write_ou_trait_evolution_summary_table(args.summary_out, report)
        if args.excluded_taxa_out:
            write_ou_trait_evolution_exclusion_table(args.excluded_taxa_out, report)
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
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "alpha": report.alpha,
                    "theta": report.theta,
                    "sigma_squared": report.sigma_squared,
                    "log_likelihood": report.log_likelihood,
                    "aic": report.aic,
                    "aicc": report.aicc,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command == "early-burst":
        report = summarize_early_burst_trait_evolution(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
        )
        if args.summary_out:
            write_early_burst_trait_evolution_summary_table(args.summary_out, report)
        if args.excluded_taxa_out:
            write_early_burst_trait_evolution_exclusion_table(
                args.excluded_taxa_out,
                report,
            )
        if args.comparison_out:
            write_early_burst_trait_evolution_comparison_table(
                args.comparison_out,
                report,
            )
        if args.profile_out:
            write_early_burst_rate_change_profile_table(args.profile_out, report)
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
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "rate_change": report.rate_change,
                    "root_state": report.root_state,
                    "sigma_squared": report.sigma_squared,
                    "log_likelihood": report.log_likelihood,
                    "aic": report.aic,
                    "aicc": report.aicc,
                    "better_model": report.better_model,
                    "identifiability_warning_count": len(
                        report.identifiability_warnings
                    ),
                    "profile_row_count": len(report.profile_rows),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command == "rate-through-time":
        report = summarize_trait_rate_through_time(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            interval_count=args.interval_count,
        )
        if args.summary_out:
            write_trait_rate_through_time_summary_table(args.summary_out, report)
        if args.intervals_out:
            write_trait_rate_through_time_interval_table(args.intervals_out, report)
        if args.excluded_taxa_out:
            write_trait_rate_through_time_exclusion_table(
                args.excluded_taxa_out,
                report,
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
                warnings=report.warnings,
                metrics={
                    "tree_taxon_count": report.tree_taxon_count,
                    "analyzed_taxon_count": report.analyzed_taxon_count,
                    "excluded_taxon_count": len(report.excluded_taxa),
                    "interval_count": report.interval_count,
                    "nonempty_interval_count": report.nonempty_interval_count,
                    "tree_depth": report.tree_depth,
                    "trend_direction": report.trend_direction,
                    "earliest_interval_rate": report.earliest_interval_rate,
                    "latest_interval_rate": report.latest_interval_rate,
                    "latest_to_earliest_rate_ratio": (
                        report.latest_to_earliest_rate_ratio
                    ),
                    "weighted_rate_slope": report.weighted_rate_slope,
                    "normalized_rate_slope": report.normalized_rate_slope,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
