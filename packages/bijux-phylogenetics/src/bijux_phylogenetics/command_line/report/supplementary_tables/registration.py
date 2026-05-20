from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports import (
    write_supplementary_batch_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .ancestral_states import (
    add_ancestral_state_supplementary_table_commands,
    run_ancestral_state_supplementary_table_command,
)
from .comparative_models import (
    add_comparative_model_supplementary_table_commands,
    run_comparative_model_supplementary_table_command,
)
from .diversification import (
    add_diversification_supplementary_table_commands,
    run_diversification_supplementary_table_command,
)
from .model_selection import (
    add_model_selection_supplementary_table_commands,
    run_model_selection_supplementary_table_command,
)
from .study_inputs import (
    add_study_input_supplementary_table_commands,
    run_study_input_supplementary_table_command,
)
from .tree_review import (
    add_tree_review_supplementary_table_commands,
    run_tree_review_supplementary_table_command,
)
def add_supplementary_table_report_commands(report_subparsers: Any) -> None:
    add_study_input_supplementary_table_commands(report_subparsers)
    add_tree_review_supplementary_table_commands(report_subparsers)
    add_model_selection_supplementary_table_commands(report_subparsers)
    add_comparative_model_supplementary_table_commands(report_subparsers)
    add_ancestral_state_supplementary_table_commands(report_subparsers)
    add_diversification_supplementary_table_commands(report_subparsers)

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


def run_supplementary_table_report_command(args: Any) -> int | None:
    study_input_result = run_study_input_supplementary_table_command(args)
    if study_input_result is not None:
        return study_input_result

    tree_review_result = run_tree_review_supplementary_table_command(args)
    if tree_review_result is not None:
        return tree_review_result

    model_selection_result = run_model_selection_supplementary_table_command(args)
    if model_selection_result is not None:
        return model_selection_result

    comparative_model_result = run_comparative_model_supplementary_table_command(args)
    if comparative_model_result is not None:
        return comparative_model_result

    ancestral_state_result = run_ancestral_state_supplementary_table_command(args)
    if ancestral_state_result is not None:
        return ancestral_state_result

    diversification_result = run_diversification_supplementary_table_command(args)
    if diversification_result is not None:
        return diversification_result

    if args.report_command == "supplementary-batch-summary-table":
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

    return None
