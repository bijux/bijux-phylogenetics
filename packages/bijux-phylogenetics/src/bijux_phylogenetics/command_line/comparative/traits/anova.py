from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.regression import (
    summarize_phylogenetic_anova,
    write_phylogenetic_anova_exclusion_table,
    write_phylogenetic_anova_group_table,
    write_phylogenetic_anova_pairwise_table,
    write_phylogenetic_anova_simulation_table,
    write_phylogenetic_anova_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_phylogenetic_anova_review_command(comparative_subparsers: Any) -> None:
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


def run_phylogenetic_anova_review_command(args: Any) -> int | None:
    if args.comparative_command != "phylogenetic-anova":
        return None

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
