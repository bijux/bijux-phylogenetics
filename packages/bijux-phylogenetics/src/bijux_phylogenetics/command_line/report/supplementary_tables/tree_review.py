from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports import (
    write_supplementary_clade_support_table,
    write_supplementary_tree_diagnostics_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_tree_review_supplementary_table_commands(report_subparsers: Any) -> None:
    report_supplementary_tree_table = report_subparsers.add_parser(
        "supplementary-tree-table",
        help="Write a supplementary tree diagnostics table with topology, support, and warning summaries.",
    )
    report_supplementary_tree_table.add_argument("--tree", required=True, type=Path)
    report_supplementary_tree_table.add_argument("--out", required=True, type=Path)
    report_supplementary_tree_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_tree_table)

    report_supplementary_clade_support_table = report_subparsers.add_parser(
        "supplementary-clade-support-table",
        help="Write a supplementary clade-support table from one reference tree and optional tree-set frequencies.",
    )
    report_supplementary_clade_support_table.add_argument(
        "--tree", required=True, type=Path
    )
    report_supplementary_clade_support_table.add_argument(
        "--comparison-tree-set", type=Path
    )
    report_supplementary_clade_support_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_clade_support_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_clade_support_table)


def run_tree_review_supplementary_table_command(args: Any) -> int | None:
    if args.report_command == "supplementary-tree-table":
        result = write_supplementary_tree_diagnostics_table(
            args.out,
            tree_path=args.tree,
        )
        inputs = [args.tree]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path],
        )
        if args.json:
            row = result.rows[0]
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=row.warnings,
                    metrics={
                        "row_count": result.row_count,
                        "tip_count": row.tip_count,
                        "supported_branch_count": row.supported_branch_count,
                        "polytomy_count": row.polytomy_count,
                        "warning_count": row.warning_count,
                        "ultrametric": row.ultrametric,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "supplementary-clade-support-table":
        result = write_supplementary_clade_support_table(
            args.out,
            tree_path=args.tree,
            comparison_tree_set_path=args.comparison_tree_set,
        )
        inputs = [args.tree]
        if args.comparison_tree_set is not None:
            inputs.append(args.comparison_tree_set)
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
                        "supported_clade_count": result.supported_clade_count,
                        "frequency_scored_clade_count": result.frequency_scored_clade_count,
                        "frequency_partial_support_count": result.frequency_partial_support_count,
                        "frequency_absent_clade_count": result.frequency_absent_clade_count,
                        "frequency_unscored_clade_count": result.frequency_unscored_clade_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    return None
