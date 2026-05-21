from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.datasets.study_inputs import inspect_metadata_table
from bijux_phylogenetics.runtime.results import build_command_result


def add_metadata_commands(subparsers: Any) -> None:
    metadata = subparsers.add_parser(
        get_command_spec("metadata").name, help=get_command_spec("metadata").summary
    )
    metadata_subparsers = metadata.add_subparsers(
        dest="metadata_command", required=True
    )
    metadata_inspect = metadata_subparsers.add_parser(
        "inspect", help="Inspect a metadata table keyed by taxon."
    )
    metadata_inspect.add_argument("table", type=Path)
    metadata_inspect.add_argument("--taxon-column")
    metadata_inspect.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(metadata_inspect)


def run_metadata_command(args: Any) -> int:
    report = inspect_metadata_table(args.table, taxon_column=args.taxon_column)
    outputs = _finalize_outputs(args, command="metadata", inputs=[args.table])
    _print_result(
        build_command_result(
            command="metadata",
            inputs=[args.table],
            outputs=outputs,
            metrics={
                "row_count": report.row_count,
                "column_count": report.column_count,
                "taxon_count": len(report.taxa),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
