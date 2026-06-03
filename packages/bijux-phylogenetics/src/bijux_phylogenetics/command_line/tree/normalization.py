from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.taxa import normalize_tree_taxa, write_taxon_mapping
from bijux_phylogenetics.runtime.results import build_command_result


def add_tree_normalization_commands(subparsers: Any) -> None:
    normalize = subparsers.add_parser(
        get_command_spec("normalize").name,
        help=get_command_spec("normalize").summary,
    )
    normalize.add_argument("tree", type=Path)
    normalize.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    normalize.add_argument("--out", required=True, type=Path)
    normalize.add_argument(
        "--json", action="store_true", help="Emit the normalization result as JSON."
    )
    _add_manifest_argument(normalize)

    normalize_taxa = subparsers.add_parser(
        get_command_spec("normalize-taxa").name,
        help=get_command_spec("normalize-taxa").summary,
    )
    normalize_taxa.add_argument("tree", type=Path)
    normalize_taxa.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    normalize_taxa.add_argument(
        "--policy", choices=("spaces-to-underscores",), required=True
    )
    normalize_taxa.add_argument("--out", required=True, type=Path)
    normalize_taxa.add_argument("--mapping-out", type=Path)
    normalize_taxa.add_argument(
        "--json", action="store_true", help="Emit the normalization result as JSON."
    )
    _add_manifest_argument(normalize_taxa)


def run_normalize_command(args: Any) -> int:
    tree = load_tree(args.tree, source_format=args.format)
    output_path = write_newick(args.out, tree)
    outputs = _finalize_outputs(
        args,
        command="normalize",
        inputs=[args.tree],
        outputs=[output_path],
    )
    if args.json:
        _print_result(
            build_command_result(
                command="normalize",
                inputs=[args.tree],
                outputs=outputs,
                metrics={"tip_count": tree.tip_count},
                data={
                    "source_format": tree.source_format,
                    "output_format": "newick",
                },
            ),
            json_output=True,
        )
    else:
        print(output_path)
    return 0


def run_normalize_taxa_command(args: Any) -> int:
    tree = load_tree(args.tree, source_format=args.format)
    normalized_tree, report = normalize_tree_taxa(tree, policy=args.policy)
    output_path = write_newick(args.out, normalized_tree)
    mapping_path = args.mapping_out or args.out.with_suffix(
        f"{args.out.suffix}.mapping.tsv"
    )
    write_taxon_mapping(mapping_path, report.renamed_taxa)
    outputs = _finalize_outputs(
        args,
        command="normalize-taxa",
        inputs=[args.tree],
        outputs=[output_path, mapping_path],
    )
    if args.json:
        _print_result(
            build_command_result(
                command="normalize-taxa",
                inputs=[args.tree],
                outputs=outputs,
                metrics={"renamed_taxa": len(report.renamed_taxa)},
                data=report,
            ),
            json_output=True,
        )
    else:
        print(output_path)
    return 0
