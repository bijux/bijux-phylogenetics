from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

from bijux_phylogenetics import __version__
from bijux_phylogenetics.command_line.registry import COMMAND_SPECS, get_command_spec
from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.diagnostics.validation import inspect_tree_path, validate_tree_path
from bijux_phylogenetics.evidence.bundles import bundle_directory
from bijux_phylogenetics.reports.service import annotate_tree_against_table, render_phylogenetics_report


def _json_ready(value: Any) -> Any:
    if is_dataclass(value):
        return {key: _json_ready(item) for key, item in asdict(value).items()}
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: _json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    return value


def _print_result(result: Any, *, json_output: bool) -> None:
    if json_output:
        print(json.dumps(_json_ready(result), indent=2, sort_keys=True))
        return
    if isinstance(result, str):
        print(result)
        return
    print(json.dumps(_json_ready(result), indent=2, sort_keys=True))


def _print_commands(*, output_format: str) -> None:
    payload = [_json_ready(spec) for spec in COMMAND_SPECS]
    if output_format == "json":
        print(json.dumps(payload, indent=2, sort_keys=True))
        return
    for command in payload:
        print(f"{command['name']}: {command['domain']} - {command['summary']}")


def build_parser() -> argparse.ArgumentParser:
    """Build the repository CLI parser."""
    parser = argparse.ArgumentParser(prog="bijux-phylogenetics")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    commands = subparsers.add_parser("commands", help="List the registered command taxonomy.")
    commands.add_argument("--format", choices=("text", "json"), default="text")

    validate = subparsers.add_parser(get_command_spec("validate").name, help=get_command_spec("validate").summary)
    validate.add_argument("tree", type=Path)
    validate.add_argument("--json", action="store_true", help="Emit the report as JSON.")

    inspect = subparsers.add_parser(get_command_spec("inspect").name, help=get_command_spec("inspect").summary)
    inspect.add_argument("tree", type=Path)
    inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")

    normalize = subparsers.add_parser(get_command_spec("normalize").name, help=get_command_spec("normalize").summary)
    normalize.add_argument("tree", type=Path)
    normalize.add_argument("--out", required=True, type=Path)

    compare = subparsers.add_parser(get_command_spec("compare").name, help=get_command_spec("compare").summary)
    compare.add_argument("left", type=Path)
    compare.add_argument("right", type=Path)
    compare.add_argument("--json", action="store_true", help="Emit the report as JSON.")

    annotate = subparsers.add_parser(get_command_spec("annotate").name, help=get_command_spec("annotate").summary)
    annotate.add_argument("tree", type=Path)
    annotate.add_argument("--metadata", required=True, type=Path)
    annotate.add_argument("--json", action="store_true", help="Emit the linkage report as JSON.")

    diagnose = subparsers.add_parser(get_command_spec("diagnose").name, help=get_command_spec("diagnose").summary)
    diagnose.add_argument("tree", type=Path)
    diagnose.add_argument("--json", action="store_true", help="Emit the report as JSON.")

    render = subparsers.add_parser(get_command_spec("render").name, help=get_command_spec("render").summary)
    render.add_argument("tree", type=Path)
    render.add_argument("--metadata", type=Path)
    render.add_argument("--out", required=True, type=Path)

    evidence = subparsers.add_parser(get_command_spec("evidence").name, help=get_command_spec("evidence").summary)
    evidence.add_argument("run_root", type=Path)
    evidence.add_argument("--out", required=True, type=Path)
    evidence.add_argument("--json", action="store_true", help="Emit the bundle report as JSON.")

    report = subparsers.add_parser(get_command_spec("report").name, help=get_command_spec("report").summary)
    report.add_argument("--tree", required=True, type=Path)
    report.add_argument("--alignment", type=Path)
    report.add_argument("--traits", type=Path)
    report.add_argument("--metadata", type=Path)
    report.add_argument("--out", required=True, type=Path)

    adapter = subparsers.add_parser(get_command_spec("adapter").name, help=get_command_spec("adapter").summary)
    adapter.add_argument("adapter_name")
    adapter.add_argument("--json", action="store_true", help="Emit the adapter report as JSON.")

    return parser


def run_command(args: Any, *, parser: argparse.ArgumentParser) -> int:
    """Run the selected command."""
    try:
        if args.command == "commands":
            _print_commands(output_format=args.format)
            return 0
        if args.command == "validate":
            _print_result(validate_tree_path(args.tree), json_output=args.json)
            return 0
        if args.command == "inspect":
            _print_result(inspect_tree_path(args.tree), json_output=args.json)
            return 0
        if args.command == "normalize":
            parser.exit(status=2, message="normalize is not implemented yet\n")
        if args.command == "diagnose":
            parser.exit(status=2, message="diagnose is not implemented yet\n")
        if args.command == "compare":
            _print_result(compare_tree_paths(args.left, args.right), json_output=args.json)
            return 0
        if args.command == "annotate":
            _print_result(annotate_tree_against_table(args.tree, args.metadata), json_output=args.json)
            return 0
        if args.command == "render":
            result = render_phylogenetics_report(tree_path=args.tree, metadata_path=args.metadata, out_path=args.out)
            print(result.output_path)
            return 0
        if args.command == "evidence":
            _print_result(bundle_directory(args.run_root, args.out), json_output=args.json)
            return 0
        if args.command == "report":
            result = render_phylogenetics_report(
                tree_path=args.tree,
                alignment_path=args.alignment,
                traits_path=args.traits,
                metadata_path=args.metadata,
                out_path=args.out,
            )
            print(result.output_path)
            return 0
        if args.command == "adapter":
            parser.exit(status=2, message=f"adapter is not implemented yet for {args.adapter_name}\n")
    except FileNotFoundError as error:
        parser.exit(status=2, message=f"{error}\n")
    except ValueError as error:
        parser.exit(status=2, message=f"{error}\n")
    except NotImplementedError as error:
        parser.exit(status=2, message=f"{error}\n")
    except Exception as error:  # pragma: no cover - defensive CLI guard
        parser.exit(status=1, message=f"unexpected error: {error}\n")

    parser.print_help(sys.stderr)
    return 2
