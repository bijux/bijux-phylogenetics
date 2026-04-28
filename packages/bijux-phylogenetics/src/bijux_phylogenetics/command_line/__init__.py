from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

from bijux_phylogenetics import __version__
from bijux_phylogenetics.command_line.registry import COMMAND_SPECS, get_command_spec
from bijux_phylogenetics.core.environment import inspect_environment
from bijux_phylogenetics.core.manifest import build_run_manifest, write_run_manifest
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.diagnostics.root_to_tip import (
    compute_root_to_tip_distances,
    diagnose_ultrametricity,
    write_root_to_tip_tsv,
)
from bijux_phylogenetics.io.fasta import build_alignment_quality_report, link_alignment_to_tree, summarise_fasta
from bijux_phylogenetics.io.fasta import (
    detect_identical_duplicate_sequences,
    detect_invalid_alignment_characters,
    detect_near_duplicate_sequences,
)
from bijux_phylogenetics.core.metadata import inspect_metadata_table
from bijux_phylogenetics.core.pruning import (
    drop_tree_taxa,
    prune_tree_to_requested_taxa,
    prune_tree_to_taxa,
    write_pruned_taxa,
)
from bijux_phylogenetics.core.traits import link_tree_to_traits, prune_traits_to_tree, validate_traits_table
from bijux_phylogenetics.core.traits import detect_missing_trait_values
from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_clade_sets,
    compare_support_values,
    compare_tree_paths,
    detect_clade_changes,
    prune_trees_to_shared_taxa,
    write_tree_comparison_table,
)
from bijux_phylogenetics.compare.reports import build_tree_comparison_report
from bijux_phylogenetics.core.demo import run_capability_demo
from bijux_phylogenetics.diagnostics.validation import diagnose_tree_path, inspect_tree_path, validate_tree_path
from bijux_phylogenetics.evidence.bundles import bundle_directory, validate_bundle
from bijux_phylogenetics.errors import EngineUnavailableError, EvidenceContractError, MetadataJoinError, PhylogeneticsError
from bijux_phylogenetics.core.taxonomy import normalize_tree_taxa, write_taxon_mapping
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.svg import render_tree_svg
from bijux_phylogenetics.reports.service import (
    annotate_tree_against_table,
    render_dataset_report,
    render_phylo_inputs_report,
    render_tree_report,
    write_annotation_report,
)
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
    if args.command == "env":
        return []
    if args.command == "metadata":
        return [args.table]
    if args.command == "traits":
        if args.traits_command == "validate":
            return [args.table]
        if args.traits_command == "missing":
            return [args.table]
        if args.traits_command == "prune":
            return [args.tree, args.table, args.out]
        return [args.tree, args.table]
    if args.command == "prune":
        inputs = [args.tree, args.out]
        if getattr(args, "keep_from", None) is not None:
            inputs.append(args.keep_from)
        if args.pruned_taxa_out is not None:
            inputs.append(args.pruned_taxa_out)
        return inputs
    if args.command == "alignment":
        if args.alignment_command == "inspect":
            return [args.alignment]
        if args.alignment_command == "composition":
            return [args.alignment]
        if args.alignment_command == "duplicates":
            return [args.alignment]
        if args.alignment_command == "invalid":
            return [args.alignment]
        if args.alignment_command == "quality":
            return [args.alignment]
        return [args.tree, args.alignment]
    if args.command in {"validate", "inspect"}:
        return [args.tree]
    if args.command == "diagnose":
        if getattr(args, "tree", None) is not None:
            return [args.tree]
        return [Path(args.target)]
    if args.command == "normalize":
        return [args.tree, args.out]
    if args.command == "normalize-taxa":
        inputs = [args.tree, args.out]
        if args.mapping_out is not None:
            inputs.append(args.mapping_out)
        return inputs
    if args.command == "compare":
        if getattr(args, "third", None) is not None:
            return [Path(args.right), Path(args.third)]
        return [Path(args.left), Path(args.right)]
    if args.command == "annotate":
        inputs = [args.tree, args.metadata]
        if args.out is not None:
            inputs.append(args.out)
        if getattr(args, "joined_out", None) is not None:
            inputs.append(args.joined_out)
        return inputs
    if args.command == "render":
        inputs = [args.tree, args.out]
        if args.metadata is not None:
            inputs.append(args.metadata)
        return inputs
    if args.command == "evidence":
        if args.evidence_command == "bundle":
            return [*args.inputs, *args.outputs]
        return [args.bundle_root]
    if args.command == "demo":
        return []
    if args.command == "report":
        if args.report_command == "tree":
            return [args.tree, args.out]
        if args.report_command == "dataset":
            inputs = [args.tree, args.metadata, args.out]
            if args.traits is not None:
                inputs.append(args.traits)
            return inputs
        if args.report_command == "phylo-inputs":
            return [args.tree, args.alignment, args.out]
        inputs = [args.tree, args.out]
        if getattr(args, "alignment", None) is not None:
            inputs.append(args.alignment)
        if getattr(args, "traits", None) is not None:
            inputs.append(args.traits)
        if getattr(args, "metadata", None) is not None:
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

    env = subparsers.add_parser(get_command_spec("env").name, help=get_command_spec("env").summary)
    env_subparsers = env.add_subparsers(dest="env_command", required=True)
    env_inspect = env_subparsers.add_parser("inspect", help="Inspect runtime dependency availability.")
    env_inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(env_inspect)

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
    traits_missing = traits_subparsers.add_parser("missing", help="List missing trait values by taxon and column.")
    traits_missing.add_argument("table", type=Path)
    traits_missing.add_argument("--taxon-column")
    traits_missing.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(traits_missing)
    traits_link = traits_subparsers.add_parser("link", help="Link tree tips to a traits table.")
    traits_link.add_argument("tree", type=Path)
    traits_link.add_argument("table", type=Path)
    traits_link.add_argument("--taxon-column")
    traits_link.add_argument("--strict", action="store_true")
    traits_link.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(traits_link)
    traits_prune = traits_subparsers.add_parser("prune", help="Prune a traits table to tree taxa.")
    traits_prune.add_argument("tree", type=Path)
    traits_prune.add_argument("table", type=Path)
    traits_prune.add_argument("--taxon-column")
    traits_prune.add_argument("--out", required=True, type=Path)
    traits_prune.add_argument("--json", action="store_true", help="Emit the pruning report as JSON.")
    _add_manifest_argument(traits_prune)

    prune = subparsers.add_parser(get_command_spec("prune").name, help=get_command_spec("prune").summary)
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

    alignment = subparsers.add_parser(get_command_spec("alignment").name, help=get_command_spec("alignment").summary)
    alignment_subparsers = alignment.add_subparsers(dest="alignment_command", required=True)
    alignment_inspect = alignment_subparsers.add_parser("inspect", help="Inspect an aligned FASTA file.")
    alignment_inspect.add_argument("alignment", type=Path)
    alignment_inspect.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_inspect)
    alignment_quality = alignment_subparsers.add_parser("quality", help="Generate a higher-level alignment quality report.")
    alignment_quality.add_argument("alignment", type=Path)
    alignment_quality.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_quality)
    alignment_composition = alignment_subparsers.add_parser(
        "composition",
        help="Inspect inferred alphabet, composition, and GC content.",
    )
    alignment_composition.add_argument("alignment", type=Path)
    alignment_composition.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_composition)
    alignment_invalid = alignment_subparsers.add_parser(
        "invalid",
        help="List alignment characters invalid for a declared alphabet.",
    )
    alignment_invalid.add_argument("alignment", type=Path)
    alignment_invalid.add_argument("--alphabet", choices=("dna", "rna", "protein"), required=True)
    alignment_invalid.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_invalid)
    alignment_duplicates = alignment_subparsers.add_parser(
        "duplicates",
        help="Report identical and near-duplicate aligned sequences.",
    )
    alignment_duplicates.add_argument("alignment", type=Path)
    alignment_duplicates.add_argument("--identity-threshold", type=float, default=0.95)
    alignment_duplicates.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(alignment_duplicates)
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
    validate.add_argument("--require-rooted", action="store_true")
    validate.add_argument("--require-ultrametric", action="store_true")
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
    compare.add_argument("left")
    compare.add_argument("right")
    compare.add_argument("third", nargs="?")
    compare.add_argument("--out", type=Path)
    compare.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(compare)

    annotate = subparsers.add_parser(get_command_spec("annotate").name, help=get_command_spec("annotate").summary)
    annotate.add_argument("tree", type=Path)
    annotate.add_argument("--metadata", required=True, type=Path)
    annotate.add_argument("--taxon-column")
    annotate.add_argument("--out", type=Path)
    annotate.add_argument("--joined-out", type=Path)
    annotate.add_argument("--json", action="store_true", help="Emit the linkage report as JSON.")
    _add_manifest_argument(annotate)

    diagnose = subparsers.add_parser(get_command_spec("diagnose").name, help=get_command_spec("diagnose").summary)
    diagnose.add_argument("target")
    diagnose.add_argument("tree", nargs="?", type=Path)
    diagnose.add_argument("--out", type=Path)
    diagnose.add_argument("--tolerance", type=float, default=1e-6)
    diagnose.add_argument("--json", action="store_true", help="Emit the report as JSON.")
    _add_manifest_argument(diagnose)

    render = subparsers.add_parser(get_command_spec("render").name, help=get_command_spec("render").summary)
    render.add_argument("tree", type=Path)
    render.add_argument("--metadata", type=Path)
    render.add_argument("--taxon-column")
    render.add_argument("--label-column")
    render.add_argument("--out", required=True, type=Path)
    render.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(render)

    evidence = subparsers.add_parser(get_command_spec("evidence").name, help=get_command_spec("evidence").summary)
    evidence_subparsers = evidence.add_subparsers(dest="evidence_command", required=True)
    evidence_bundle = evidence_subparsers.add_parser(
        "bundle",
        help="Bundle explicit phylogenetics inputs and outputs as evidence.",
    )
    evidence_bundle.add_argument("--inputs", nargs="+", required=True, type=Path)
    evidence_bundle.add_argument("--outputs", nargs="+", required=True, type=Path)
    evidence_bundle.add_argument("--out", required=True, type=Path)
    evidence_bundle.add_argument("--json", action="store_true", help="Emit the bundle report as JSON.")
    _add_manifest_argument(evidence_bundle)
    evidence_validate = evidence_subparsers.add_parser("validate", help="Validate an existing evidence bundle.")
    evidence_validate.add_argument("bundle_root", type=Path)
    evidence_validate.add_argument("--json", action="store_true", help="Emit the validation report as JSON.")
    _add_manifest_argument(evidence_validate)

    report = subparsers.add_parser(get_command_spec("report").name, help=get_command_spec("report").summary)
    report_subparsers = report.add_subparsers(dest="report_command", required=True)
    report_tree = report_subparsers.add_parser("tree", help="Render a deterministic single-tree HTML report.")
    report_tree.add_argument("tree", type=Path)
    report_tree.add_argument("--out", required=True, type=Path)
    report_tree.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(report_tree)
    report_dataset = report_subparsers.add_parser("dataset", help="Render a tree plus table dataset HTML report.")
    report_dataset.add_argument("--tree", required=True, type=Path)
    report_dataset.add_argument("--metadata", required=True, type=Path)
    report_dataset.add_argument("--traits", type=Path)
    report_dataset.add_argument("--out", required=True, type=Path)
    report_dataset.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(report_dataset)
    report_phylo_inputs = report_subparsers.add_parser(
        "phylo-inputs",
        help="Render a tree plus alignment HTML input report.",
    )
    report_phylo_inputs.add_argument("--tree", required=True, type=Path)
    report_phylo_inputs.add_argument("--alignment", required=True, type=Path)
    report_phylo_inputs.add_argument("--out", required=True, type=Path)
    report_phylo_inputs.add_argument("--json", action="store_true", help="Emit the report build result as JSON.")
    _add_manifest_argument(report_phylo_inputs)

    demo = subparsers.add_parser(get_command_spec("demo").name, help=get_command_spec("demo").summary)
    demo_subparsers = demo.add_subparsers(dest="demo_command", required=True)
    demo_run = demo_subparsers.add_parser("run", help="Run the repository capability demo workflow.")
    demo_run.add_argument("--out", required=True, type=Path)
    demo_run.add_argument("--json", action="store_true", help="Emit the demo result as JSON.")
    _add_manifest_argument(demo_run)

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
        if args.command == "env":
            report = inspect_environment()
            outputs = _finalize_outputs(args, command="env", inputs=[])
            _print_result(
                build_command_result(
                    command="env",
                    inputs=[],
                    outputs=outputs,
                    metrics={"dependency_count": len(report.dependencies)},
                    data=report,
                ),
                json_output=args.json,
            )
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
        if args.command == "validate":
            report = validate_tree_path(
                args.tree,
                source_format=args.format,
                allow_duplicates=args.allow_duplicates,
                strict=args.strict,
                allow_negative_branch_lengths=args.allow_negative_branches,
                require_rooted=args.require_rooted,
                require_ultrametric=args.require_ultrametric,
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
            if args.keep_from is not None:
                tree, report = prune_tree_to_taxa(args.tree, args.keep_from, taxon_column=args.taxon_column)
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
            if args.alignment_command == "quality":
                report = build_alignment_quality_report(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        warnings=report.warnings,
                        metrics={
                            "invalid_character_count": len(report.invalid_characters),
                            "composition_outlier_count": len(report.composition_outliers),
                            "duplicate_group_count": len(report.duplicate_sequence_groups),
                            "near_duplicate_count": len(report.near_duplicate_pairs),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "composition":
                report = summarise_fasta(args.alignment)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "alphabet": report.inferred_alphabet,
                            "sequence_count": report.sequence_count,
                            "gc_sequence_count": len(report.per_sequence_gc_content),
                        },
                        data={
                            "path": report.path,
                            "inferred_alphabet": report.inferred_alphabet,
                            "nucleotide_composition": report.nucleotide_composition,
                            "amino_acid_composition": report.amino_acid_composition,
                            "per_sequence_gc_content": report.per_sequence_gc_content,
                            "whole_alignment_gc_content": report.whole_alignment_gc_content,
                        },
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "invalid":
                report = detect_invalid_alignment_characters(args.alignment, alphabet=args.alphabet)
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={"invalid_character_count": len(report), "alphabet": args.alphabet},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.alignment_command == "duplicates":
                duplicates = detect_identical_duplicate_sequences(args.alignment)
                near_duplicates = detect_near_duplicate_sequences(
                    args.alignment,
                    identity_threshold=args.identity_threshold,
                )
                outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
                _print_result(
                    build_command_result(
                        command="alignment",
                        inputs=[args.alignment],
                        outputs=outputs,
                        metrics={
                            "duplicate_group_count": len(duplicates),
                            "near_duplicate_count": len(near_duplicates),
                            "identity_threshold": args.identity_threshold,
                        },
                        data={
                            "duplicate_sequence_groups": duplicates,
                            "near_duplicate_pairs": near_duplicates,
                        },
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
                    warnings=report.warnings,
                    metrics={
                        "tip_count": report.tip_count,
                        "node_count": report.node_count,
                        "internal_node_count": report.internal_node_count,
                        "edge_count": report.edge_count,
                        "clade_count": report.clade_count,
                        "is_binary": report.is_binary,
                        "polytomy_count": report.polytomy_count,
                        "branch_length_status": report.branch_length_status,
                        "is_ultrametric": report.is_ultrametric,
                        "tree_diameter": report.tree_diameter,
                        "colless_imbalance_index": report.colless_imbalance_index,
                        "sackin_imbalance_index": report.sackin_imbalance_index,
                        "tree_quality_score": report.tree_quality_score,
                        "zero_length_branch_count": report.zero_length_branch_count,
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
            if args.target == "distances":
                if args.tree is None:
                    parser.exit(status=2, message="diagnose distances requires a tree path\n")
                report = compute_root_to_tip_distances(args.tree)
                outputs: list[Path | str] = []
                if args.out is not None:
                    output_path = write_root_to_tip_tsv(args.out, report)
                    outputs.append(output_path)
                outputs = _finalize_outputs(args, command="diagnose", inputs=[args.tree], outputs=outputs)
                _print_result(
                    build_command_result(
                        command="diagnose",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={"tip_count": len(report.distances)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.target == "ultrametric":
                if args.tree is None:
                    parser.exit(status=2, message="diagnose ultrametric requires a tree path\n")
                report = diagnose_ultrametricity(args.tree, tolerance=args.tolerance)
                outputs = _finalize_outputs(args, command="diagnose", inputs=[args.tree])
                _print_result(
                    build_command_result(
                        command="diagnose",
                        inputs=[args.tree],
                        outputs=outputs,
                        metrics={
                            "tip_count": report.tip_count,
                            "tolerance": report.tolerance,
                            "max_deviation": report.max_deviation,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            tree_path = Path(args.target) if args.tree is None else args.tree
            report = diagnose_tree_path(tree_path)
            outputs = _finalize_outputs(args, command="diagnose", inputs=[tree_path])
            _print_result(
                build_command_result(
                    command="diagnose",
                    inputs=[tree_path],
                    outputs=outputs,
                    warnings=report.validation.warnings,
                    metrics={
                        "tip_count": report.inspection.tip_count,
                        "polytomy_count": report.validation.polytomy_count,
                        "cherry_count": report.inspection.cherry_count,
                        "tree_diameter": report.inspection.tree_diameter,
                        "tree_quality_score": report.inspection.tree_quality_score,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "compare":
            if args.left == "report":
                if args.third is None:
                    parser.exit(status=2, message="compare report requires two tree paths\n")
                if args.out is None:
                    parser.exit(status=2, message="compare report requires --out\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                report = build_tree_comparison_report(left_path, right_path, out_path=args.out)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path], outputs=[args.out])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={"shared_taxa": len(report.topology.shared_taxa)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "support":
                if args.third is None:
                    parser.exit(status=2, message="compare support requires two tree paths\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                report = compare_support_values(left_path, right_path)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={"shared_clades": len(report.shared_clades)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "clades":
                if args.third is None:
                    parser.exit(status=2, message="compare clades requires two tree paths\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                report = compare_clade_sets(left_path, right_path)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={
                            "shared_clades": len(report.shared_clades),
                            "left_only_clades": len(report.left_only_clades),
                            "right_only_clades": len(report.right_only_clades),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "prune":
                if args.third is None:
                    parser.exit(status=2, message="compare prune requires two tree paths\n")
                if args.out is None:
                    parser.exit(status=2, message="compare prune requires --out as an output directory\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                pruned_left, pruned_right, report = prune_trees_to_shared_taxa(left_path, right_path)
                args.out.mkdir(parents=True, exist_ok=True)
                left_out = write_newick(args.out / "left-shared.nwk", pruned_left)
                right_out = write_newick(args.out / "right-shared.nwk", pruned_right)
                outputs = _finalize_outputs(
                    args,
                    command="compare",
                    inputs=[left_path, right_path],
                    outputs=[left_out, right_out],
                )
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={"shared_taxa": len(report.shared_taxa)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "changes":
                if args.third is None:
                    parser.exit(status=2, message="compare changes requires two tree paths\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                report = detect_clade_changes(left_path, right_path)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={
                            "lost_clades": len(report.lost_clades),
                            "gained_clades": len(report.gained_clades),
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "branch-lengths":
                if args.third is None:
                    parser.exit(status=2, message="compare branch-lengths requires two tree paths\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                report = compare_branch_lengths(left_path, right_path)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={"shared_splits": len(report.shared_splits)},
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            if args.left == "table":
                if args.third is None:
                    parser.exit(status=2, message="compare table requires two tree paths\n")
                if args.out is None:
                    parser.exit(status=2, message="compare table requires --out\n")
                left_path = Path(args.right)
                right_path = Path(args.third)
                output_path = write_tree_comparison_table(args.out, left_path, right_path)
                outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path], outputs=[output_path])
                _print_result(
                    build_command_result(
                        command="compare",
                        inputs=[left_path, right_path],
                        outputs=outputs,
                        metrics={"table_rows": sum(1 for _ in output_path.read_text(encoding='utf-8').splitlines()[1:])},
                        data={"table_path": output_path},
                    ),
                    json_output=args.json,
                )
                return 0
            left_path = Path(args.left)
            right_path = Path(args.right)
            report = compare_tree_paths(left_path, right_path)
            outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
            _print_result(
                build_command_result(
                    command="compare",
                    inputs=[left_path, right_path],
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
            report = annotate_tree_against_table(args.tree, args.metadata, taxon_column=args.taxon_column)
            outputs: list[Path | str] = []
            if args.out is not None:
                outputs.append(write_annotation_report(args.out, report))
            if args.joined_out is not None:
                table = load_taxon_table(args.metadata, taxon_column=args.taxon_column)
                outputs.append(
                    write_taxon_rows(
                        args.joined_out,
                        columns=["taxon", "matched", *[column for column in table.columns if column != table.taxon_column]],
                        rows=[
                            {
                                "taxon": row.taxon,
                                "matched": str(row.matched).lower(),
                                **{column: row.values.get(column, "") for column in table.columns if column != table.taxon_column},
                            }
                            for row in report.joined_rows
                        ],
                    )
                )
            outputs = _finalize_outputs(args, command="annotate", inputs=[args.tree, args.metadata], outputs=outputs)
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
            labels: dict[str, str] | None = None
            if args.metadata is not None and args.label_column is not None:
                table = load_taxon_table(args.metadata, taxon_column=args.taxon_column)
                if args.label_column not in table.columns:
                    raise MetadataJoinError(f"metadata table does not contain label column '{args.label_column}'")
                labels = {row[table.taxon_column]: row[args.label_column] for row in table.rows if row[args.label_column]}
            result = render_tree_svg(args.tree, out_path=args.out, labels=labels)
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
                        warnings=result.missing_metadata_labels,
                        metrics={"tip_count": result.tip_count, "label_count": result.label_count},
                        data=result,
                    ),
                    json_output=True,
                )
                return 0
            print(result.output_path)
            return 0
        if args.command == "evidence":
            if args.evidence_command == "bundle":
                report = bundle_directory(args.inputs, args.outputs, args.out)
                inputs = [*args.inputs, *args.outputs]
                outputs = _finalize_outputs(args, command="evidence", inputs=inputs, outputs=[args.out])
                _print_result(
                    build_command_result(
                        command="evidence",
                        inputs=inputs,
                        outputs=outputs,
                        metrics={
                            "file_count": report.file_count,
                            "input_file_count": report.input_file_count,
                            "output_file_count": report.output_file_count,
                        },
                        data=report,
                    ),
                    json_output=args.json,
                )
                return 0
            report = validate_bundle(args.bundle_root)
            if not report.valid:
                raise EvidenceContractError(f"evidence bundle validation failed with {len(report.mismatches)} mismatch(es)")
            outputs = _finalize_outputs(args, command="evidence", inputs=[args.bundle_root])
            _print_result(
                build_command_result(
                    command="evidence",
                    inputs=[args.bundle_root],
                    outputs=outputs,
                    metrics={"file_count": report.file_count},
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        if args.command == "demo":
            if args.demo_command == "run":
                result = run_capability_demo(args.out)
                outputs = _finalize_outputs(
                    args,
                    command="demo",
                    inputs=[],
                    outputs=[
                        result.tree_report,
                        result.dataset_report,
                        result.phylo_inputs_report,
                        result.comparison_report,
                        result.evidence_bundle,
                        result.capability_summary,
                    ],
                )
                if args.json:
                    _print_result(
                        build_command_result(
                            command="demo",
                            inputs=[],
                            outputs=outputs,
                            metrics={"artifact_count": 6},
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_root)
                return 0
            raise NotImplementedError(f"unsupported demo command: {args.demo_command}")
        if args.command == "report":
            if args.report_command == "tree":
                result = render_tree_report(tree_path=args.tree, out_path=args.out)
                outputs = _finalize_outputs(args, command="report", inputs=[args.tree], outputs=[result.output_path])
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=[args.tree],
                            outputs=outputs,
                            warnings=result.validation.warnings + result.inspection.warnings,
                            metrics={"tip_count": result.inspection.tip_count},
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            if args.report_command == "dataset":
                result = render_dataset_report(
                    tree_path=args.tree,
                    metadata_path=args.metadata,
                    traits_path=args.traits,
                    out_path=args.out,
                )
                inputs = [args.tree, args.metadata]
                if args.traits is not None:
                    inputs.append(args.traits)
                outputs = _finalize_outputs(args, command="report", inputs=inputs, outputs=[result.output_path])
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=inputs,
                            outputs=outputs,
                            warnings=result.validation.warnings + result.inspection.warnings,
                            metrics={"tip_count": result.inspection.tip_count, "linked_taxa": result.metadata_linkage.linked_taxa},
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            if args.report_command == "phylo-inputs":
                result = render_phylo_inputs_report(
                    tree_path=args.tree,
                    alignment_path=args.alignment,
                    out_path=args.out,
                )
                inputs = [args.tree, args.alignment]
                outputs = _finalize_outputs(args, command="report", inputs=inputs, outputs=[result.output_path])
                if args.json:
                    _print_result(
                        build_command_result(
                            command="report",
                            inputs=inputs,
                            outputs=outputs,
                            warnings=result.validation.warnings + result.inspection.warnings,
                            metrics={
                                "tip_count": result.inspection.tip_count,
                                "alignment_length": result.alignment.alignment_length,
                                "linked_taxa": result.alignment_linkage.linked_taxa,
                            },
                            data=result,
                        ),
                        json_output=True,
                    )
                    return 0
                print(result.output_path)
                return 0
            raise NotImplementedError(f"unsupported report command: {args.report_command}")
        if args.command == "adapter":
            raise EngineUnavailableError(f"adapter is not available in this runtime: {args.adapter_name}")
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
