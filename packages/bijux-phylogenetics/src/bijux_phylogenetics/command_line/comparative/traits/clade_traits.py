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
from bijux_phylogenetics.runtime.results import build_command_result


def add_clade_trait_review_command(comparative_subparsers: Any) -> None:
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


def run_clade_trait_review_command(args: Any) -> int | None:
    if args.comparative_command != "clade-traits":
        return None

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
