from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.datasets.study_inputs import (
    check_tree_and_trait_taxon_names,
    detect_missing_trait_values,
    link_tree_to_traits,
    load_taxon_table,
    prune_traits_to_tree,
    validate_traits_table,
    write_taxon_rows,
    write_tree_trait_name_mismatch_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_traits_commands(subparsers: Any) -> None:
    traits = subparsers.add_parser(
        get_command_spec("traits").name, help=get_command_spec("traits").summary
    )
    traits_subparsers = traits.add_subparsers(dest="traits_command", required=True)
    traits_validate = traits_subparsers.add_parser(
        "validate", help="Validate a traits table keyed by taxon."
    )
    traits_validate.add_argument("table", type=Path)
    traits_validate.add_argument("--taxon-column")
    traits_validate.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(traits_validate)

    traits_missing = traits_subparsers.add_parser(
        "missing", help="List missing trait values by taxon and column."
    )
    traits_missing.add_argument("table", type=Path)
    traits_missing.add_argument("--taxon-column")
    traits_missing.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(traits_missing)

    traits_link = traits_subparsers.add_parser(
        "link", help="Link tree tips to a traits table."
    )
    traits_link.add_argument("tree", type=Path)
    traits_link.add_argument("table", type=Path)
    traits_link.add_argument("--taxon-column")
    traits_link.add_argument("--strict", action="store_true")
    traits_link.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(traits_link)

    traits_name_check = traits_subparsers.add_parser(
        "name-check", help="Report tree-versus-trait taxon mismatches."
    )
    traits_name_check.add_argument("tree", type=Path)
    traits_name_check.add_argument("table", type=Path)
    traits_name_check.add_argument("--taxon-column")
    traits_name_check.add_argument("--out", type=Path)
    traits_name_check.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(traits_name_check)

    traits_prune = traits_subparsers.add_parser(
        "prune", help="Prune a traits table to tree taxa."
    )
    traits_prune.add_argument("tree", type=Path)
    traits_prune.add_argument("table", type=Path)
    traits_prune.add_argument("--taxon-column")
    traits_prune.add_argument("--out", required=True, type=Path)
    traits_prune.add_argument(
        "--json", action="store_true", help="Emit the pruning report as JSON."
    )
    _add_manifest_argument(traits_prune)


def run_traits_command(args: Any) -> int:
    if args.traits_command == "validate":
        report = validate_traits_table(args.table, taxon_column=args.taxon_column)
        outputs = _finalize_outputs(args, command="traits", inputs=[args.table])
        _print_result(
            build_command_result(
                command="traits",
                inputs=[args.table],
                outputs=outputs,
                metrics={
                    "row_count": report.row_count,
                    "trait_column_count": len(report.trait_columns),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.traits_command == "missing":
        report = detect_missing_trait_values(args.table, taxon_column=args.taxon_column)
        outputs = _finalize_outputs(args, command="traits", inputs=[args.table])
        _print_result(
            build_command_result(
                command="traits",
                inputs=[args.table],
                outputs=outputs,
                metrics={"missing_value_count": len(report.missing_values)},
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.traits_command == "name-check":
        report = check_tree_and_trait_taxon_names(
            args.tree,
            args.table,
            taxon_column=args.taxon_column,
        )
        declared_outputs: list[Path] = []
        if args.out is not None:
            declared_outputs.append(
                write_tree_trait_name_mismatch_table(args.out, report)
            )
        outputs = _finalize_outputs(
            args,
            command="traits",
            inputs=[args.tree, args.table],
            outputs=declared_outputs,
        )
        _print_result(
            build_command_result(
                command="traits",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "tree_not_data_count": len(report.tree_not_data),
                    "data_not_tree_count": len(report.data_not_tree),
                    "compatible": report.compatible,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.traits_command == "prune":
        rows, report = prune_traits_to_tree(
            args.tree,
            args.table,
            taxon_column=args.taxon_column,
        )
        table = load_taxon_table(args.table, taxon_column=args.taxon_column)
        output_path = write_taxon_rows(args.out, columns=table.columns, rows=rows)
        outputs = _finalize_outputs(
            args,
            command="traits",
            inputs=[args.tree, args.table],
            outputs=[output_path],
        )
        _print_result(
            build_command_result(
                command="traits",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "original_row_count": report.original_row_count,
                    "kept_taxa": len(report.kept_taxa),
                    "removed_taxa": len(report.removed_taxa),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    report = link_tree_to_traits(
        args.tree,
        args.table,
        taxon_column=args.taxon_column,
        strict=args.strict,
    )
    outputs = _finalize_outputs(args, command="traits", inputs=[args.tree, args.table])
    _print_result(
        build_command_result(
            command="traits",
            inputs=[args.tree, args.table],
            outputs=outputs,
            metrics={
                "tree_taxa": report.tree_taxa,
                "trait_taxa": report.trait_taxa,
                "linked_taxa": report.linked_taxa,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
