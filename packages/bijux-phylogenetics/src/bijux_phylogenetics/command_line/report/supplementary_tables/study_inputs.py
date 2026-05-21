from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports import (
    write_supplementary_alignment_diagnostics_table,
    write_supplementary_taxon_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_study_input_supplementary_table_commands(report_subparsers: Any) -> None:
    report_supplementary_taxon_table = report_subparsers.add_parser(
        "supplementary-taxon-table",
        help="Write a supplementary taxon table with IDs, metadata, traits, and exclusion evidence.",
    )
    report_supplementary_taxon_table.add_argument("--tree", required=True, type=Path)
    report_supplementary_taxon_table.add_argument(
        "--metadata", required=True, type=Path
    )
    report_supplementary_taxon_table.add_argument("--traits", required=True, type=Path)
    report_supplementary_taxon_table.add_argument("--alignment", type=Path)
    report_supplementary_taxon_table.add_argument("--filtered-alignment", type=Path)
    report_supplementary_taxon_table.add_argument("--inference-tree", type=Path)
    report_supplementary_taxon_table.add_argument("--reported-taxa", type=Path)
    report_supplementary_taxon_table.add_argument("--tip-dates", type=Path)
    report_supplementary_taxon_table.add_argument("--calibrations", type=Path)
    report_supplementary_taxon_table.add_argument("--out", required=True, type=Path)
    report_supplementary_taxon_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_taxon_table)

    report_supplementary_alignment_table = report_subparsers.add_parser(
        "supplementary-alignment-table",
        help="Write a supplementary alignment diagnostics table with optional filtering status.",
    )
    report_supplementary_alignment_table.add_argument(
        "--alignment", required=True, type=Path
    )
    report_supplementary_alignment_table.add_argument("--filtered-alignment", type=Path)
    report_supplementary_alignment_table.add_argument("--out", required=True, type=Path)
    report_supplementary_alignment_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_alignment_table)


def run_study_input_supplementary_table_command(args: Any) -> int | None:
    if args.report_command == "supplementary-taxon-table":
        result = write_supplementary_taxon_table(
            args.out,
            tree_path=args.tree,
            metadata_path=args.metadata,
            traits_path=args.traits,
            alignment_path=args.alignment,
            filtered_alignment_path=args.filtered_alignment,
            inference_tree_path=args.inference_tree,
            reported_taxa_path=args.reported_taxa,
            tip_dates_path=args.tip_dates,
            calibration_path=args.calibrations,
        )
        inputs = [args.tree, args.metadata, args.traits]
        if args.alignment is not None:
            inputs.append(args.alignment)
        if args.filtered_alignment is not None:
            inputs.append(args.filtered_alignment)
        if args.inference_tree is not None:
            inputs.append(args.inference_tree)
        if args.reported_taxa is not None:
            inputs.append(args.reported_taxa)
        if args.tip_dates is not None:
            inputs.append(args.tip_dates)
        if args.calibrations is not None:
            inputs.append(args.calibrations)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "analysis_included_count": result.analysis_included_count,
                        "analysis_excluded_count": result.analysis_excluded_count,
                        "reporting_retained_count": result.reporting_retained_count,
                        "reporting_dropped_count": result.reporting_dropped_count,
                        "metadata_column_count": result.metadata_column_count,
                        "trait_column_count": result.trait_column_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-alignment-table":
        result = write_supplementary_alignment_diagnostics_table(
            args.out,
            alignment_path=args.alignment,
            filtered_alignment_path=args.filtered_alignment,
        )
        inputs = [args.alignment]
        if args.filtered_alignment is not None:
            inputs.append(args.filtered_alignment)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=[],
                    metrics={
                        "row_count": result.row_count,
                        "retained_sequence_count": result.retained_sequence_count,
                        "removed_sequence_count": result.removed_sequence_count,
                        "filtered_only_sequence_count": result.filtered_only_sequence_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    return None
