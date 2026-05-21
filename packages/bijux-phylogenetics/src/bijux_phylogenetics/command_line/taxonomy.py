from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _parse_labelled_run,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.taxa import (
    audit_tree_taxon_synonyms,
    build_taxon_audit_report,
    build_taxon_stability_report,
    build_taxon_workflow_loss_report,
    export_tree_accepted_names,
    inspect_tree_taxon_namespaces,
    inspect_tree_taxon_rank_consistency,
    load_taxon_run_source,
    resolve_tree_taxon_synonyms,
    write_accepted_name_mapping,
    write_synonym_resolution_mapping,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_taxonomy_commands(subparsers: Any) -> None:
    taxonomy = subparsers.add_parser(
        get_command_spec("taxonomy").name, help=get_command_spec("taxonomy").summary
    )
    taxonomy_subparsers = taxonomy.add_subparsers(
        dest="taxonomy_command", required=True
    )

    taxonomy_synonyms = taxonomy_subparsers.add_parser(
        "synonyms",
        help="Audit a tree against a configurable taxon synonym table.",
    )
    taxonomy_synonyms.add_argument("tree", type=Path)
    taxonomy_synonyms.add_argument("--synonym-table", required=True, type=Path)
    taxonomy_synonyms.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    taxonomy_synonyms.add_argument(
        "--json", action="store_true", help="Emit the audit as JSON."
    )
    _add_manifest_argument(taxonomy_synonyms)

    taxonomy_resolve = taxonomy_subparsers.add_parser(
        "resolve-synonyms",
        help="Resolve tree tip synonyms to accepted labels with a reversible mapping artifact.",
    )
    taxonomy_resolve.add_argument("tree", type=Path)
    taxonomy_resolve.add_argument("--synonym-table", required=True, type=Path)
    taxonomy_resolve.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    taxonomy_resolve.add_argument(
        "--resolution-policy", choices=("reject-ambiguous",), default="reject-ambiguous"
    )
    taxonomy_resolve.add_argument("--out", required=True, type=Path)
    taxonomy_resolve.add_argument("--mapping-out", type=Path)
    taxonomy_resolve.add_argument(
        "--json", action="store_true", help="Emit the resolution report as JSON."
    )
    _add_manifest_argument(taxonomy_resolve)

    taxonomy_namespaces = taxonomy_subparsers.add_parser(
        "namespaces",
        help="Classify tree tip labels into accession, species, sample, isolate, or user-defined namespaces.",
    )
    taxonomy_namespaces.add_argument("tree", type=Path)
    taxonomy_namespaces.add_argument(
        "--format", choices=("newick", "nexus", "phyloxml")
    )
    taxonomy_namespaces.add_argument(
        "--json", action="store_true", help="Emit the namespace report as JSON."
    )
    _add_manifest_argument(taxonomy_namespaces)

    taxonomy_ranks = taxonomy_subparsers.add_parser(
        "rank-consistency",
        help="Audit whether tree labels mix species, genus, sample, accession, or population naming levels.",
    )
    taxonomy_ranks.add_argument("tree", type=Path)
    taxonomy_ranks.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    taxonomy_ranks.add_argument(
        "--json", action="store_true", help="Emit the rank audit as JSON."
    )
    _add_manifest_argument(taxonomy_ranks)

    taxonomy_accepted = taxonomy_subparsers.add_parser(
        "accepted-names",
        help="Export raw tree labels to accepted names through a configurable synonym table.",
    )
    taxonomy_accepted.add_argument("tree", type=Path)
    taxonomy_accepted.add_argument("--synonym-table", required=True, type=Path)
    taxonomy_accepted.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    taxonomy_accepted.add_argument("--out", type=Path)
    taxonomy_accepted.add_argument(
        "--json", action="store_true", help="Emit the accepted-name export as JSON."
    )
    _add_manifest_argument(taxonomy_accepted)

    taxonomy_audit = taxonomy_subparsers.add_parser(
        "audit",
        help="Build a reviewer-readable taxon audit across namespaces, ranks, synonyms, and mapping conflicts.",
    )
    taxonomy_audit.add_argument("tree", type=Path)
    taxonomy_audit.add_argument("--synonym-table", type=Path)
    taxonomy_audit.add_argument("--format", choices=("newick", "nexus", "phyloxml"))
    taxonomy_audit.add_argument(
        "--json", action="store_true", help="Emit the taxon audit as JSON."
    )
    _add_manifest_argument(taxonomy_audit)

    taxonomy_loss = taxonomy_subparsers.add_parser(
        "loss",
        help="Trace where taxa were lost across tree, alignment, metadata, traits, inference, and reporting stages.",
    )
    taxonomy_loss.add_argument("tree", type=Path)
    taxonomy_loss.add_argument("--metadata", required=True, type=Path)
    taxonomy_loss.add_argument("--traits", required=True, type=Path)
    taxonomy_loss.add_argument("--alignment", type=Path)
    taxonomy_loss.add_argument("--filtered-alignment", type=Path)
    taxonomy_loss.add_argument("--inference-tree", type=Path)
    taxonomy_loss.add_argument("--reported-taxa", type=Path)
    taxonomy_loss.add_argument(
        "--json", action="store_true", help="Emit the workflow loss report as JSON."
    )
    _add_manifest_argument(taxonomy_loss)

    taxonomy_stability = taxonomy_subparsers.add_parser(
        "stability",
        help="Compare taxon retention across multiple named workflow artifacts.",
    )
    taxonomy_stability.add_argument(
        "--run",
        action="append",
        required=True,
        metavar="LABEL=PATH",
        help="Named run source in the form label=/path/to/tree_or_alignment_or_table",
    )
    taxonomy_stability.add_argument(
        "--json", action="store_true", help="Emit the stability report as JSON."
    )
    _add_manifest_argument(taxonomy_stability)


def run_taxonomy_command(args: Any) -> int:
    if args.taxonomy_command == "synonyms":
        tree = load_tree(args.tree, source_format=args.format)
        report = audit_tree_taxon_synonyms(tree, args.synonym_table)
        outputs = _finalize_outputs(
            args, command="taxonomy", inputs=[args.tree, args.synonym_table]
        )
        _print_result(
            build_command_result(
                command="taxonomy",
                inputs=[args.tree, args.synonym_table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "candidate_count": len(report.candidates),
                    "ambiguous_mapping_count": len(report.ambiguous_mappings),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.taxonomy_command == "resolve-synonyms":
        tree = load_tree(args.tree, source_format=args.format)
        resolved_tree, report = resolve_tree_taxon_synonyms(
            tree,
            synonym_table_path=args.synonym_table,
            resolution_policy=args.resolution_policy,
        )
        output_path = write_newick(args.out, resolved_tree)
        mapping_path = args.mapping_out or args.out.with_suffix(
            f"{args.out.suffix}.synonyms.tsv"
        )
        write_synonym_resolution_mapping(mapping_path, report)
        outputs = _finalize_outputs(
            args,
            command="taxonomy",
            inputs=[args.tree, args.synonym_table],
            outputs=[output_path, mapping_path],
        )
        _print_result(
            build_command_result(
                command="taxonomy",
                inputs=[args.tree, args.synonym_table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "renamed_taxa": len(report.renamed_taxa),
                    "ambiguous_mapping_count": len(report.ambiguous_mappings),
                    "duplicate_resolved_label_count": len(
                        report.duplicate_resolved_labels
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.taxonomy_command == "namespaces":
        tree = load_tree(args.tree, source_format=args.format)
        report = inspect_tree_taxon_namespaces(tree)
        outputs = _finalize_outputs(args, command="taxonomy", inputs=[args.tree])
        _print_result(
            build_command_result(
                command="taxonomy",
                inputs=[args.tree],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "namespace_count": len(report.namespace_counts),
                    "mixed_namespaces": report.mixed_namespaces,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.taxonomy_command == "rank-consistency":
        tree = load_tree(args.tree, source_format=args.format)
        report = inspect_tree_taxon_rank_consistency(tree)
        outputs = _finalize_outputs(args, command="taxonomy", inputs=[args.tree])
        _print_result(
            build_command_result(
                command="taxonomy",
                inputs=[args.tree],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "rank_count": len(report.rank_counts),
                    "mixed_ranks": report.mixed_ranks,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.taxonomy_command == "accepted-names":
        tree = load_tree(args.tree, source_format=args.format)
        report = export_tree_accepted_names(tree, args.synonym_table)
        output_paths: list[Path] = []
        if args.out is not None:
            output_paths.append(write_accepted_name_mapping(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="taxonomy",
            inputs=[args.tree, args.synonym_table],
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="taxonomy",
                inputs=[args.tree, args.synonym_table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "row_count": len(report.rows),
                    "ambiguous_count": sum(
                        1 for row in report.rows if row.status == "ambiguous"
                    ),
                    "resolved_count": sum(
                        1 for row in report.rows if row.status == "resolved"
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        if not args.json and args.out is not None:
            print(args.out)
        return 0

    if args.taxonomy_command == "audit":
        tree = load_tree(args.tree, source_format=args.format)
        report = build_taxon_audit_report(tree, synonym_table_path=args.synonym_table)
        inputs = [args.tree]
        if args.synonym_table is not None:
            inputs.append(args.synonym_table)
        outputs = _finalize_outputs(args, command="taxonomy", inputs=inputs)
        _print_result(
            build_command_result(
                command="taxonomy",
                inputs=inputs,
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "tree_tip_count": report.tree_tip_count,
                    "status": report.status,
                    "mapping_conflict_count": len(report.mapping_conflicts.rows),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.taxonomy_command == "loss":
        report = build_taxon_workflow_loss_report(
            args.tree,
            args.metadata,
            args.traits,
            alignment_path=args.alignment,
            filtered_alignment_path=args.filtered_alignment,
            inference_tree_path=args.inference_tree,
            reported_taxa_path=args.reported_taxa,
        )
        inputs = [args.tree, args.metadata, args.traits]
        for optional_path in (
            args.alignment,
            args.filtered_alignment,
            args.inference_tree,
            args.reported_taxa,
        ):
            if optional_path is not None:
                inputs.append(optional_path)
        outputs = _finalize_outputs(args, command="taxonomy", inputs=inputs)
        _print_result(
            build_command_result(
                command="taxonomy",
                inputs=inputs,
                outputs=outputs,
                metrics={
                    "taxon_count": len(report.rows),
                    "loss_stage_count": len(report.loss_stage_counts),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    run_sources = [
        load_taxon_run_source(label=label, path=path)
        for label, path in (_parse_labelled_run(raw) for raw in args.run)
    ]
    report = build_taxon_stability_report(run_sources)
    inputs = [source.path for source in run_sources]
    outputs = _finalize_outputs(args, command="taxonomy", inputs=inputs)
    _print_result(
        build_command_result(
            command="taxonomy",
            inputs=inputs,
            outputs=outputs,
            metrics={
                "source_count": len(report.sources),
                "stable_taxa": len(report.stable_taxa),
                "unstable_taxa": len(report.unstable_taxa),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
