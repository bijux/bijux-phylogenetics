from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.engines import render_inference_workflow_report
from bijux_phylogenetics.runtime.results import build_command_result


def add_adapter_reporting_commands(adapter_subparsers: Any) -> None:
    adapter_report = adapter_subparsers.add_parser(
        "report", help="Render an HTML report from an engine workflow manifest."
    )
    adapter_report.add_argument("manifest_path", type=Path)
    adapter_report.add_argument("--out", required=True, type=Path)
    adapter_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(adapter_report)


def run_adapter_reporting_command(args: Any) -> int | None:
    if args.adapter_command != "report":
        return None

    report = render_inference_workflow_report(
        manifest_path=args.manifest_path,
        out_path=args.out,
    )
    outputs = _finalize_outputs(
        args,
        command="adapter",
        inputs=[args.manifest_path],
        outputs=[report.output_path],
    )
    _print_result(
        build_command_result(
            command="adapter",
            inputs=[args.manifest_path],
            outputs=outputs,
            metrics={"warning_count": report.warning_count},
            data=report,
        ),
        json_output=args.json,
    )
    return 0
