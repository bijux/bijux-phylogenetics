from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports import write_supplementary_model_selection_table
from bijux_phylogenetics.runtime.results import build_command_result


def add_model_selection_supplementary_table_commands(
    report_subparsers: Any,
) -> None:
    report_supplementary_model_selection_table = report_subparsers.add_parser(
        "supplementary-model-selection-table",
        help="Write a supplementary model-selection table from parsed IQ-TREE report artifacts.",
    )
    report_supplementary_model_selection_table.add_argument(
        "--iqtree-report", required=True, type=Path
    )
    report_supplementary_model_selection_table.add_argument(
        "--model-sidecar", type=Path
    )
    report_supplementary_model_selection_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_model_selection_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_model_selection_table)


def run_model_selection_supplementary_table_command(args: Any) -> int | None:
    if args.report_command != "supplementary-model-selection-table":
        return None

    result = write_supplementary_model_selection_table(
        args.out,
        iqtree_report_path=args.iqtree_report,
        model_sidecar_path=args.model_sidecar,
    )
    inputs = [args.iqtree_report]
    if args.model_sidecar is not None:
        inputs.append(args.model_sidecar)
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
                    "candidate_count": result.candidate_count,
                    "selected_model": result.selected_model,
                    "selected_criterion": result.selected_criterion,
                },
                data=result,
            ),
            json_output=True,
        )
        return 0
    print(result.output_path)
    return 0
