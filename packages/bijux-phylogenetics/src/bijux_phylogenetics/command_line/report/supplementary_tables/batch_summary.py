from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports import write_supplementary_batch_summary_table
from bijux_phylogenetics.runtime.results import build_command_result


def add_batch_summary_supplementary_table_commands(report_subparsers: Any) -> None:
    report_supplementary_batch_summary_table = report_subparsers.add_parser(
        "supplementary-batch-summary-table",
        help="Write a supplementary batch summary table from a written workflow bundle.",
    )
    report_supplementary_batch_summary_table.add_argument(
        "--workflow-bundle-root", required=True, type=Path
    )
    report_supplementary_batch_summary_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_batch_summary_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_batch_summary_table)


def run_batch_summary_supplementary_table_command(args: Any) -> int | None:
    if args.report_command != "supplementary-batch-summary-table":
        return None

    result = write_supplementary_batch_summary_table(
        args.out,
        workflow_bundle_root=args.workflow_bundle_root,
    )
    inputs = [args.workflow_bundle_root]
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
                    "dataset_row_count": result.dataset_row_count,
                    "variant_row_count": result.variant_row_count,
                    "workflow_status": result.workflow_status,
                    "warning_count": result.warning_count,
                },
                data=result,
            ),
            json_output=True,
        )
        return 0
    print(result.output_path)
    return 0
