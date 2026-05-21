from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.reports.service import (
    annotate_tree_against_table,
    write_annotation_report,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_annotate_command(subparsers: Any) -> None:
    annotate = subparsers.add_parser(
        get_command_spec("annotate").name, help=get_command_spec("annotate").summary
    )
    annotate.add_argument("tree", type=Path)
    annotate.add_argument("--metadata", required=True, type=Path)
    annotate.add_argument("--taxon-column")
    annotate.add_argument("--out", type=Path)
    annotate.add_argument("--joined-out", type=Path)
    annotate.add_argument(
        "--json", action="store_true", help="Emit the linkage report as JSON."
    )
    _add_manifest_argument(annotate)


def run_annotate_command(args: Any) -> int:
    report = annotate_tree_against_table(
        args.tree, args.metadata, taxon_column=args.taxon_column
    )
    outputs: list[Path | str] = []
    if args.out is not None:
        outputs.append(write_annotation_report(args.out, report))
    if args.joined_out is not None:
        table = load_taxon_table(args.metadata, taxon_column=args.taxon_column)
        outputs.append(
            write_taxon_rows(
                args.joined_out,
                columns=[
                    "taxon",
                    "matched",
                    *[
                        column
                        for column in table.columns
                        if column != table.taxon_column
                    ],
                ],
                rows=[
                    {
                        "taxon": row.taxon,
                        "matched": str(row.matched).lower(),
                        **{
                            column: row.values.get(column, "")
                            for column in table.columns
                            if column != table.taxon_column
                        },
                    }
                    for row in report.joined_rows
                ],
            )
        )
    outputs = _finalize_outputs(
        args,
        command="annotate",
        inputs=[args.tree, args.metadata],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="annotate",
            inputs=[args.tree, args.metadata],
            outputs=outputs,
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
