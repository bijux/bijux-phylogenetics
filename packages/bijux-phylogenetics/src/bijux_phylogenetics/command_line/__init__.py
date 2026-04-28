from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

from bijux_phylogenetics import __version__
from bijux_phylogenetics.command_line.registry import COMMAND_SPECS, get_command_spec
from bijux_phylogenetics.core.manifest import build_run_manifest, write_run_manifest
from bijux_phylogenetics.io.fasta import link_alignment_to_tree, summarise_fasta
from bijux_phylogenetics.core.metadata import inspect_metadata_table
from bijux_phylogenetics.core.pruning import prune_tree_to_taxa, write_pruned_taxa
from bijux_phylogenetics.core.traits import link_tree_to_traits, validate_traits_table
from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.diagnostics.validation import diagnose_tree_path, inspect_tree_path, validate_tree_path
from bijux_phylogenetics.evidence.bundles import bundle_directory
from bijux_phylogenetics.errors import PhylogeneticsError
from bijux_phylogenetics.core.taxonomy import normalize_tree_taxa, write_taxon_mapping
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


def _add_manifest_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--manifest", type=Path, help="Write a reproducibility manifest to this JSON path.")


def _finalize_outputs(
    args: Any,
    *,
    command: str,
    inputs: list[Path | str],
    outputs: list[Path | str] | None = None,
) -> list[Path | str]:
    finalized_outputs = list(outputs or [])
    manifest_path = getattr(args, "manifest", None)
    if manifest_path is None:
        return finalized_outputs
    manifest = build_run_manifest(
        command=command,
        arguments=list(getattr(args, "_argv", [])),
        input_paths=inputs,
        output_paths=finalized_outputs,
    )
    write_run_manifest(manifest_path, manifest)
    finalized_outputs.append(manifest_path)
    return finalized_outputs


def _command_inputs(args: Any) -> list[Path | str]:
    if args.command == "commands":
        return []
    if args.command == "metadata":
        return [args.table]
    if args.command == "traits":
        if args.traits_command == "validate":
            return [args.table]
        return [args.tree, args.table]
    if args.command == "prune":
        inputs = [args.tree, args.keep_from, args.out]
        if args.pruned_taxa_out is not None:
            inputs.append(args.pruned_taxa_out)
        return inputs
    if args.command == "alignment":
        if args.alignment_command == "inspect":
            return [args.alignment]
        return [args.tree, args.alignment]
    if args.command in {"validate", "inspect", "diagnose"}:
        return [args.tree]
    if args.command == "normalize":
        return [args.tree, args.out]
    if args.command == "normalize-taxa":
        inputs = [args.tree, args.out]
        if args.mapping_out is not None:
            inputs.append(args.mapping_out)
        return inputs
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

    metadata = subparsers.add_parser(get_command_spec("metadata").name, help=get_command_spec("metadata").summary)
    metadata_subparsers = metadata.add_subparsers(dest="metadata_command", required=True)
    metadata_inspect = metadata_subparsers.add_parser("inspect", help="Inspect a metadata table keyed by taxon.")
    metadata_inspect.add_argument("table", type=Path)
    metadata_inspect.add_argument("--taxon-column")
    metadata_inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(metadata_inspect)

    traits = subparsers.add_parser(get_command_spec("traits").name, help=get_command_spec("traits").summary)
    traits_subparsers = traits.add_subparsers(dest="traits_command", required=True)
    traits_validate = traits_subparsers.add_parser("validate", help="Validate a traits table keyed by taxon.")
    traits_validate.add_argument("table", type=Path)
    traits_validate.add_argument("--taxon-column")
    traits_validate.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(traits_validate)
    traits_link = traits_subparsers.add_parser("link", help="Link tree tips to a traits table.")
    traits_link.add_argument("tree", type=Path)
    traits_link.add_argument("table", type=Path)
    traits_link.add_argument("--taxon-column")
    traits_link.add_argument("--strict", action="store_true")
    traits_link.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(traits_link)

    prune = subparsers.add_parser(get_command_spec("prune").name, help=get_command_spec("prune").summary)
    prune.add_argument("tree", type=Path)
    prune.add_argument("--keep-from", required=True, type=Path)
    prune.add_argument("--taxon-column")
    prune.add_argument("--out", required=True, type=Path)
    prune.add_argument("--pruned-taxa-out", type=Path)
    prune.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(prune)

    alignment = subparsers.add_parser(get_command_spec("alignment").name, help=get_command_spec("alignment").summary)
    alignment_subparsers = alignment.add_subparsers(dest="alignment_command", required=True)
    alignment_inspect = alignment_subparsers.add_parser("inspect", help="Inspect an aligned FASTA file.")
    alignment_inspect.add_argument("alignment", type=Path)
    alignment_inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_inspect)
    alignment_link = alignment_subparsers.add_parser("link", help="Link tree tips to an aligned FASTA file.")
    alignment_link.add_argument("tree", type=Path)
    alignment_link.add_argument("alignment", type=Path)
    alignment_link.add_argument("--strict", action="store_true")
    alignment_link.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_link)

    validate = subparsers.add_parser(get_command_spec("validate").name, help=get_command_spec("validate").summary)
    validate.add_argument("tree", type=Path)
    validate.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    validate.add_argument("--allow-duplicates", action="store_true")
    validate.add_argument("--allow-negative-branches", action="store_true")
    validate.add_argument("--strict", action="store_true")
    validate.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(validate)

    inspect = subparsers.add_parser(get_command_spec("inspect").name, help=get_command_spec("inspect").summary)
    inspect.add_argument("tree", type=Path)
    inspect.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(inspect)

    normalize = subparsers.add_parser(get_command_spec("normalize").name, help=get_command_spec("normalize").summary)
    normalize.add_argument("tree", type=Path)
    normalize.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    normalize.add_argument("--out", required=True, type=Path)
    normalize.add_argument("--json", action="store_true", help="Emit the normalization result as JSON.")
    _add_manifest_argument(normalize)

    normalize_taxa = subparsers.add_parser(
        get_command_spec("normalize-taxa").name,
        help=get_command_spec("normalize-taxa").summary,
    )
    normalize_taxa.add_argument("tree", type=Path)
    normalize_taxa.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    normalize_taxa.add_argument("--policy", choices=("spaces-to-underscores",), required=True)
    normalize_taxa.add_argument("--out", required=True, type=Path)
    normalize_taxa.add_argument("--mapping-out", type=Path)
    normalize_taxa.add_argument("--json", action="store_true", help="Emit the normalization result as JSON.")
    _add_manifest_argument(normalize_taxa)

    compare = subparsers.add_parser(get_command_spec("compare").name, help=get_command_spec("compare").summary)
    compare.add_argument("left", type=Path)
    compare.add_argument("right", type=Path)
    compare.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(compare)

    annotate = subparsers.add_parser(get_command_spec("annotate").name, help=get_command_spec("annotate").summary)
    annotate.add_argument("tree", type=Path)
    annotate.add_argument("--metadata", required=True, type=Path)
    annotate.add_argument("--json", action="store_true", help="Emit the linkage report as JSON.")
    _add_manifest_argument(annotate)

    diagnose = subparsers.add_parser(get_command_spec("diagnose").name, help=get_command_spec("diagnose").summary)
    diagnose.add_argument("tree", type=Path)
    diagnose.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(diagnose)

    render = subparsers.add_parser(get_command_spec("render").name, help=get_command_spec("render").summary)
    render.add_argument("tree", type=Path)
    render.add_argument("--metadata", type=Path)
    render.add_argument("--out", required=True, type=Path)
    render.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(render)

    evidence = subparsers.add_parser(get_command_spec("evidence").name, help=get_command_spec("evidence").summary)
    evidence.add_argument("run_root", type=Path)
    evidence.add_argument("--out", required=True, type=Path)
    evidence.add_argument("--json", action="store_true", help="Emit the bundle report as JSON.")
    _add_manifest_argument(evidence)

    report = subparsers.add_parser(get_command_spec("report").name, help=get_command_spec("report").summary)
    report.add_argument("--tree", required=True, type=Path)
    report.add_argument("--alignment", type=Path)
    report.add_argument("--traits", type=Path)
    report.add_argument("--metadata", type=Path)
    report.add_argument("--out", required=True, type=Path)
    report.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(report)

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
        if args.command == "metadata":
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
        if args.command == "traits":
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
        if args.command == "validate":
            report = validate_tree_path(
                args.tree,
                source_format=args.format,
                allow_duplicates=args.allow_duplicates,
                strict=args.strict,
                allow_negative_branch_lengths=args.allow_negative_branches,
            )
            outputs = _finalize_outputs(args, command="validate", inputs=[args.tree])
            _print_result(
                build_command_result(
                    command="validate",
                    inputs=[args.tree],
                    outputs=outputs,
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
        if args.command == "prune":
            tree, report = prune_tree_to_taxa(args.tree, args.keep_from, taxon_column=args.taxon_column)
            output_path = write_newick(args.out, tree)
            pruned_taxa_path = args.pruned_taxa_out or args.out.with_name("pruned_taxa.tsv")
            write_pruned_taxa(pruned_taxa_path, report.removed_taxa)
            outputs = _finalize_outputs(
                args,
                command="prune",
                inputs=[args.tree, args.keep_from],
                outputs=[output_path, pruned_taxa_path],
            )
            _print_result(
                build_command_result(
                    command="prune",
                    inputs=[args.tree, args.keep_from],
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
        if args.command == "alignment":
            if args.alignment_command == "inspect":
                report = summarise_fasta(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "sequence_count": report.sequence_count,
                            "alignment_length": report.alignment_length,
                            "variable_site_count": report.variable_site_count,
                            "parsimony_informative_site_count": report.parsimony_informative_site_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            report = link_alignment_to_tree(args.tree, args.alignment, strict=args.strict)
            outputs = _finalize_outputs(args, command="alignment", inputs=[args.tree, args.alignment])
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
        if args.command == "inspect":
            report = inspect_tree_path(args.tree, source_format=args.format)
            outputs = _finalize_outputs(args, command="inspect", inputs=[args.tree])
            _print_result(
                build_command_result(
                    command="inspect",
                    inputs=[args.tree],
                    outputs=outputs,
                    metrics={
                        "tip_count": report.tip_count,
                        "node_count": report.node_count,
                        "internal_node_count": report.internal_node_count,
                        "edge_count": report.edge_count,
                        "clade_count": report.clade_count,
                        "is_binary": report.is_binary,
                        "polytomy_count": report.polytomy_count,
                        "branch_length_status": report.branch_length_status,
                        "cherry_count": report.cherry_count,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "normalize":
            tree = load_tree(args.tree, source_format=args.format)
            output_path = write_newick(args.out, tree)
            outputs = _finalize_outputs(args, command="normalize", inputs=[args.tree], outputs=[output_path])
            if args.json:
                _print_result(
                    build_command_result(
                        command="normalize",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={"tip_count": tree.tip_count},
                        data={"source_format": tree.source_format, "output_format": "newick"},
                    ),
                    json_output=True,
                )
            else:
                print(output_path)
            return 0
        if args.command == "normalize-taxa":
            tree = load_tree(args.tree, source_format=args.format)
            normalized_tree, report = normalize_tree_taxa(tree, policy=args.policy)
            output_path = write_newick(args.out, normalized_tree)
            mapping_path = args.mapping_out or args.out.with_suffix(f"{args.out.suffix}.mapping.tsv")
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
        if args.command == "diagnose":
            report = diagnose_tree_path(args.tree)
            outputs = _finalize_outputs(args, command="diagnose", inputs=[args.tree])
            _print_result(
                build_command_result(
                    command="diagnose",
                    inputs=[args.tree],
                    outputs=outputs,
                    warnings=report.validation.warnings,
                    metrics={
                        "tip_count": report.inspection.tip_count,
                        "polytomy_count": report.validation.polytomy_count,
                        "cherry_count": report.inspection.cherry_count,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "compare":
            report = compare_tree_paths(args.left, args.right)
            outputs = _finalize_outputs(args, command="compare", inputs=[args.left, args.right])
            _print_result(
                build_command_result(
                    command="compare",
                    inputs=[args.left, args.right],
                    outputs=outputs,
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
            outputs = _finalize_outputs(args, command="annotate", inputs=[args.tree, args.metadata])
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
        if args.command == "render":
            result = render_phylogenetics_report(tree_path=args.tree, metadata_path=args.metadata, out_path=args.out)
            inputs = [args.tree]
            if args.metadata is not None:
                inputs.append(args.metadata)
            outputs = _finalize_outputs(args, command="render", inputs=inputs, outputs=[result.output_path])
            if args.json:
                _print_result(
                    build_command_result(
                        command="render",
                        inputs=inputs,
                        outputs=outputs,
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
            outputs = _finalize_outputs(args, command="evidence", inputs=[args.run_root], outputs=[args.out])
            _print_result(
                build_command_result(
                    command="evidence",
                    inputs=[args.run_root],
                    outputs=outputs,
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
            inputs = [args.tree]
            if args.alignment is not None:
                inputs.append(args.alignment)
            if args.traits is not None:
                inputs.append(args.traits)
            if args.metadata is not None:
                inputs.append(args.metadata)
            outputs = _finalize_outputs(args, command="report", inputs=inputs, outputs=[result.output_path])
            if args.json:
                _print_result(
                    build_command_result(
                        command="report",
                        inputs=inputs,
                        outputs=outputs,
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
