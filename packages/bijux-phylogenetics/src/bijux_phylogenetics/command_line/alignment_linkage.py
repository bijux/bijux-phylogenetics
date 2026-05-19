from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.fasta.records import link_alignment_to_tree
from bijux_phylogenetics.runtime.results import build_command_result


def add_alignment_linkage_commands(alignment_subparsers: Any) -> None:
    alignment_link = alignment_subparsers.add_parser(
        "link",
        help="Link tree tips to an aligned FASTA file.",
    )
    alignment_link.add_argument("tree", type=Path)
    alignment_link.add_argument("alignment", type=Path)
    alignment_link.add_argument("--strict", action="store_true")
    alignment_link.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_link)


def run_alignment_linkage_command(args: Any) -> int | None:
    if args.alignment_command != "link":
        return None

    report = link_alignment_to_tree(
        args.tree,
        args.alignment,
        strict=args.strict,
    )
    outputs = _finalize_outputs(
        args,
        command="alignment",
        inputs=[args.tree, args.alignment],
    )
    _print_result(
        build_command_result(
            command="alignment",
            inputs=[args.tree, args.alignment],
            outputs=outputs,
            metrics={
                "tree_taxa": report.tree_taxa,
                "alignment_ids": report.alignment_ids,
                "linked_taxa": report.linked_taxa,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
