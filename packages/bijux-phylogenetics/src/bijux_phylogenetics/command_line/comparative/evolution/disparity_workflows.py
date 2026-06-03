from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _split_csv_values,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.disparity import (
    render_disparity_through_time_svg,
    summarize_continuous_clade_disparity,
    summarize_disparity_through_time,
    write_continuous_clade_disparity_summary_table,
    write_continuous_clade_disparity_table,
    write_disparity_through_time_bin_table,
    write_disparity_through_time_curve_table,
    write_disparity_through_time_exclusion_table,
    write_disparity_through_time_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_disparity_workflow_comparative_evolution_commands(
    comparative_subparsers: Any,
) -> None:
    comparative_dtt = comparative_subparsers.add_parser(
        "dtt",
        help="Summarize geiger-style disparity through time for one continuous trait matrix.",
    )
    comparative_dtt.add_argument("tree", type=Path)
    comparative_dtt.add_argument("table", type=Path)
    comparative_dtt.add_argument(
        "--traits",
        required=True,
        help="Comma-delimited continuous trait columns used as the disparity matrix.",
    )
    comparative_dtt.add_argument("--taxon-column")
    comparative_dtt.add_argument(
        "--time-bin-count",
        type=int,
        help="Optional equal-width time-bin count used to summarize the raw DTT curve.",
    )
    comparative_dtt.add_argument(
        "--summary-out",
        type=Path,
        help="Write one disparity-through-time summary ledger as TSV or CSV.",
    )
    comparative_dtt.add_argument(
        "--curve-out",
        type=Path,
        help="Write one raw disparity-through-time curve ledger as TSV or CSV.",
    )
    comparative_dtt.add_argument(
        "--clades-out",
        type=Path,
        help="Write one internal-clade disparity ledger as TSV or CSV.",
    )
    comparative_dtt.add_argument(
        "--bins-out",
        type=Path,
        help="Write one equal-width time-bin disparity ledger as TSV or CSV.",
    )
    comparative_dtt.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for the DTT fit as TSV or CSV.",
    )
    comparative_dtt.add_argument(
        "--svg-out",
        type=Path,
        help="Render one SVG disparity-through-time figure.",
    )
    comparative_dtt.add_argument(
        "--json",
        action="store_true",
        help="Emit the disparity-through-time report as JSON.",
    )
    _add_manifest_argument(comparative_dtt)

    comparative_disparity = comparative_subparsers.add_parser(
        "disparity",
        help="Summarize geiger-style continuous clade disparity for one rooted tree.",
    )
    comparative_disparity.add_argument("tree", type=Path)
    comparative_disparity.add_argument("table", type=Path)
    comparative_disparity.add_argument(
        "--traits",
        required=True,
        help="Comma-delimited continuous trait columns used as the disparity matrix.",
    )
    comparative_disparity.add_argument("--taxon-column")
    comparative_disparity.add_argument(
        "--summary-out",
        type=Path,
        help="Write one clade disparity summary ledger as TSV or CSV.",
    )
    comparative_disparity.add_argument(
        "--clades-out",
        type=Path,
        help="Write one internal-clade disparity ledger as TSV or CSV.",
    )
    comparative_disparity.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for clade disparity as TSV or CSV.",
    )
    comparative_disparity.add_argument(
        "--json",
        action="store_true",
        help="Emit the clade disparity report as JSON.",
    )
    _add_manifest_argument(comparative_disparity)


def run_disparity_workflow_comparative_evolution_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    del parser
    if args.comparative_command == "dtt":
        report = summarize_disparity_through_time(
            args.tree,
            args.table,
            trait_columns=_split_csv_values(args.traits),
            taxon_column=args.taxon_column,
            time_bin_count=args.time_bin_count,
        )
        rendered_point_count = 0
        if args.summary_out:
            write_disparity_through_time_summary_table(args.summary_out, report)
        if args.curve_out:
            write_disparity_through_time_curve_table(args.curve_out, report)
        if args.clades_out:
            write_continuous_clade_disparity_table(args.clades_out, report)
        if args.bins_out:
            write_disparity_through_time_bin_table(args.bins_out, report)
        if args.excluded_taxa_out:
            write_disparity_through_time_exclusion_table(
                args.excluded_taxa_out,
                report,
            )
        if args.svg_out:
            rendered_point_count = render_disparity_through_time_svg(
                args.svg_out,
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
                    "trait_column_count": len(report.trait_columns),
                    "curve_point_count": len(report.curve_rows),
                    "time_bin_count": len(report.time_bin_rows),
                    "root_disparity": report.root_disparity,
                    "relative_scaling_applied": report.relative_scaling_applied,
                    "rendered_point_count": rendered_point_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command != "disparity":
        return None
    report = summarize_continuous_clade_disparity(
        args.tree,
        args.table,
        trait_columns=_split_csv_values(args.traits),
        taxon_column=args.taxon_column,
    )
    if args.summary_out:
        write_continuous_clade_disparity_summary_table(args.summary_out, report)
    if args.clades_out:
        write_continuous_clade_disparity_table(args.clades_out, report)
    if args.excluded_taxa_out:
        write_disparity_through_time_exclusion_table(
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
                "trait_column_count": len(report.trait_columns),
                "clade_count": len(report.clade_rows),
                "root_disparity": report.root_disparity,
                "minimum_clade_disparity": report.minimum_clade_disparity,
                "maximum_clade_disparity": report.maximum_clade_disparity,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
