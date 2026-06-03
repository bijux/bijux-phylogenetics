from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _adapter_version_args,
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.engines import read_engine_version
from bijux_phylogenetics.runtime.results import build_command_result


def add_adapter_inspection_commands(adapter_subparsers: Any) -> None:
    adapter_inspect = adapter_subparsers.add_parser(
        "inspect", help="Report external engine version metadata."
    )
    adapter_inspect.add_argument(
        "engine_name", choices=("mafft", "trimal", "iqtree", "FastTree", "MrBayes")
    )
    adapter_inspect.add_argument("--executable", type=str)
    adapter_inspect.add_argument(
        "--json", action="store_true", help="Emit the adapter report as JSON."
    )
    _add_manifest_argument(adapter_inspect)


def run_adapter_inspection_command(args: Any) -> int | None:
    if args.adapter_command != "inspect":
        return None

    executable = args.executable or args.engine_name
    report = read_engine_version(
        args.engine_name,
        executable,
        version_args=_adapter_version_args(args.engine_name),
    )
    outputs = _finalize_outputs(args, command="adapter", inputs=[args.engine_name])
    _print_result(
        build_command_result(
            command="adapter",
            inputs=[args.engine_name],
            outputs=outputs,
            metrics={"version_line_count": len(report.text.splitlines())},
            data=report,
        ),
        json_output=args.json,
    )
    return 0
