from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.phylo.pruning import (
    drop_tree_taxa,
    prune_tree_to_requested_taxa,
    prune_tree_to_taxa,
    write_pruned_taxa,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_prune_command(subparsers: Any) -> None:
    prune = subparsers.add_parser(
        get_command_spec("prune").name, help=get_command_spec("prune").summary
    )
    prune.add_argument("tree", type=Path)
    prune_targets = prune.add_mutually_exclusive_group(required=True)
    prune_targets.add_argument("--keep-from", type=Path)
    prune_targets.add_argument("--taxa", nargs="+")
    prune_targets.add_argument("--exclude-taxa", nargs="+")
    prune.add_argument("--taxon-column")
    prune.add_argument("--out", required=True, type=Path)
    prune.add_argument("--pruned-taxa-out", type=Path)
    prune.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(prune)


def run_prune_command(args: Any) -> int:
    if args.keep_from is not None:
        tree, report = prune_tree_to_taxa(
            args.tree,
            args.keep_from,
            taxon_column=args.taxon_column,
        )
        prune_inputs = [args.tree, args.keep_from]
    elif args.exclude_taxa is not None:
        tree, report = drop_tree_taxa(args.tree, list(args.exclude_taxa))
        prune_inputs = [args.tree]
    else:
        tree, report = prune_tree_to_requested_taxa(args.tree, list(args.taxa))
        prune_inputs = [args.tree]

    output_path = write_newick(args.out, tree)
    pruned_taxa_path = args.pruned_taxa_out or args.out.with_name("pruned_taxa.tsv")
    write_pruned_taxa(pruned_taxa_path, report.removed_taxa)
    outputs = _finalize_outputs(
        args,
        command="prune",
        inputs=prune_inputs,
        outputs=[output_path, pruned_taxa_path],
    )
    _print_result(
        build_command_result(
            command="prune",
            inputs=prune_inputs,
            outputs=outputs,
            metrics={
                "kept_taxa": len(report.kept_taxa),
                "removed_taxa": len(report.removed_taxa),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
