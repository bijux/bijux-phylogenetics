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
from bijux_phylogenetics.errors import PhylogeneticsError
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.reports.service import annotate_tree_against_table, render_phylogenetics_report
from bijux_phylogenetics.results import build_command_result, build_error_result


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
    payload = build_command_result(
        command="commands",
        inputs=[],
        outputs=[],
        metrics={"command_count": len(COMMAND_SPECS)},
        data={"commands": list(COMMAND_SPECS)},
    )
    if output_format == "json":
        print(json.dumps(_json_ready(payload), indent=2, sort_keys=True))
        return
    for command in _json_ready(payload.data)["commands"]:
        print(f"{command['name']}: {command['domain']} - {command['summary']}")


def _json_requested(args: Any) -> bool:
    return bool(getattr(args, "json", False) or getattr(args, "format", "") == "json")


def _command_inputs(args: Any) -> list[Path | str]:
    if args.command == "commands":
        return []
    if args.command in {"validate", "inspect", "diagnose"}:
        return [args.tree]
    if args.command == "normalize":
        return [args.tree, args.out]
    if args.command == "compare":
        return [args.left, args.right]
    if args.command == "annotate":
        return [args.tree, args.metadata]
    if args.command == "render":
        inputs = [args.tree, args.out]
        if args.metadata is not None:
            inputs.append(args.metadata)
        return inputs
    if args.command == "evidence":
        return [args.run_root, args.out]
    if args.command == "report":
        inputs = [args.tree, args.out]
        if args.alignment is not None:
            inputs.append(args.alignment)
        if args.traits is not None:
            inputs.append(args.traits)
        if args.metadata is not None:
            inputs.append(args.metadata)
        return inputs
    if args.command == "adapter":
        return [args.adapter_name]
    return []


def build_parser() -> argparse.ArgumentParser:
    """Build the repository CLI parser."""
    parser = argparse.ArgumentParser(prog="bijux-phylogenetics")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    commands = subparsers.add_parser("commands", help="List the registered command taxonomy.")
    commands.add_argument("--format", choices=("text", "json"), default="text")

    validate = subparsers.add_parser(get_command_spec("validate").name, help=get_command_spec("validate").summary)
    validate.add_argument("tree", type=Path)
    validate.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    validate.add_argument("--json", action="store_true", help="Emit the report as JSON.")

    inspect = subparsers.add_parser(get_command_spec("inspect").name, help=get_command_spec("inspect").summary)
    inspect.add_argument("tree", type=Path)
    inspect.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")

    normalize = subparsers.add_parser(get_command_spec("normalize").name, help=get_command_spec("normalize").summary)
    normalize.add_argument("tree", type=Path)
    normalize.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    normalize.add_argument("--out", required=True, type=Path)
    normalize.add_argument("--json", action="store_true", help="Emit the normalization result as JSON.")

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
    render.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")

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
    report.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")

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
            report = validate_tree_path(args.tree, source_format=args.format)
            _print_result(
                build_command_result(
                    command="validate",
                    inputs=[args.tree],
                    warnings=report.warnings,
                    metrics={
                        "tip_count": report.tip_count,
                        "internal_node_count": report.internal_node_count,
                        "polytomy_count": report.polytomy_count,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "inspect":
            report = inspect_tree_path(args.tree, source_format=args.format)
            _print_result(
                build_command_result(
                    command="inspect",
                    inputs=[args.tree],
                    metrics={
                        "tip_count": report.tip_count,
                        "internal_node_count": report.internal_node_count,
                        "edge_count": report.edge_count,
                        "is_binary": report.is_binary,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "normalize":
            tree = load_tree(args.tree, source_format=args.format)
            output_path = write_newick(args.out, tree)
            if args.json:
                _print_result(
                    build_command_result(
                        command="normalize",
                        inputs=[args.tree],
                        outputs=[output_path],
                        metrics={"tip_count": tree.tip_count},
                        data={"source_format": tree.source_format, "output_format": "newick"},
                    ),
                    json_output=True,
                )
            else:
                print(output_path)
            return 0
        if args.command == "diagnose":
            parser.exit(status=2, message="diagnose is not implemented yet\n")
        if args.command == "compare":
            report = compare_tree_paths(args.left, args.right)
            _print_result(
                build_command_result(
                    command="compare",
                    inputs=[args.left, args.right],
                    metrics={
                        "shared_taxa": len(report.shared_taxa),
                        "robinson_foulds_distance": report.robinson_foulds_distance,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "annotate":
            report = annotate_tree_against_table(args.tree, args.metadata)
            _print_result(
                build_command_result(
                    command="annotate",
                    inputs=[args.tree, args.metadata],
                    metrics={
                        "tree_taxa": report.tree_taxa,
                        "table_rows": report.table_rows,
                        "linked_taxa": report.linked_taxa,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "render":
            result = render_phylogenetics_report(tree_path=args.tree, metadata_path=args.metadata, out_path=args.out)
            if args.json:
                inputs = [args.tree]
                if args.metadata is not None:
                    inputs.append(args.metadata)
                _print_result(
                    build_command_result(
                        command="render",
                        inputs=inputs,
                        outputs=[result.output_path],
                        warnings=result.validation.warnings,
                        metrics={"tip_count": result.inspection.tip_count},
                        data=result,
                    ),
                    json_output=True,
                )
                return 0
            print(result.output_path)
            return 0
        if args.command == "evidence":
            report = bundle_directory(args.run_root, args.out)
            _print_result(
                build_command_result(
                    command="evidence",
                    inputs=[args.run_root],
                    outputs=[args.out],
                    metrics={"file_count": report.file_count},
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "report":
            result = render_phylogenetics_report(
                tree_path=args.tree,
                alignment_path=args.alignment,
                traits_path=args.traits,
                metadata_path=args.metadata,
                out_path=args.out,
            )
            if args.json:
                inputs = [args.tree]
                if args.alignment is not None:
                    inputs.append(args.alignment)
                if args.traits is not None:
                    inputs.append(args.traits)
                if args.metadata is not None:
                    inputs.append(args.metadata)
                _print_result(
                    build_command_result(
                        command="report",
                        inputs=inputs,
                        outputs=[result.output_path],
                        warnings=result.validation.warnings,
                        metrics={"tip_count": result.inspection.tip_count},
                        data=result,
                    ),
                    json_output=True,
                )
                return 0
            print(result.output_path)
            return 0
        if args.command == "adapter":
            parser.exit(status=2, message=f"adapter is not implemented yet for {args.adapter_name}\n")
    except PhylogeneticsError as error:
        if _json_requested(args):
            _print_result(
                build_error_result(command=args.command, inputs=_command_inputs(args), error=error),
                json_output=True,
            )
            return 2
        parser.exit(status=2, message=f"{error.code}: {error.message}\n")
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
