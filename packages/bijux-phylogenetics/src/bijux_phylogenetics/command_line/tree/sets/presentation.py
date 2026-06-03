from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports.service import render_tree_uncertainty_report
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.trees import (
    build_tree_set_uncertainty_method_report,
    write_bootstrap_tree_set_artifacts,
    write_tree_set_uncertainty_methods_summary_text,
)
from bijux_phylogenetics.trees.uncertainty import (
    build_tree_set_uncertainty_figure_package,
)


def add_tree_set_presentation_commands(tree_set_subparsers: Any) -> None:
    tree_set_bootstrap_summary = tree_set_subparsers.add_parser(
        "bootstrap-summary",
        help="Summarize bootstrap replicate trees into consensus and instability artifacts.",
    )
    tree_set_bootstrap_summary.add_argument("tree_set", type=Path)
    tree_set_bootstrap_summary.add_argument("--out-dir", required=True, type=Path)
    tree_set_bootstrap_summary.add_argument(
        "--prefix", default="bootstrap-tree-set", help="Prefix for written artifacts."
    )
    tree_set_bootstrap_summary.add_argument(
        "--consensus-threshold", type=float, default=0.5
    )
    tree_set_bootstrap_summary.add_argument(
        "--robust-support-threshold", type=float, default=0.9
    )
    tree_set_bootstrap_summary.add_argument(
        "--json", action="store_true", help="Emit the bootstrap summary report as JSON."
    )
    _add_manifest_argument(tree_set_bootstrap_summary)

    tree_set_methods_summary = tree_set_subparsers.add_parser(
        "methods-summary",
        help="Write reviewer-facing Markdown methods text for one tree-set uncertainty analysis.",
    )
    tree_set_methods_summary.add_argument("tree_set", type=Path)
    tree_set_methods_summary.add_argument("--out", required=True, type=Path)
    tree_set_methods_summary.add_argument(
        "--json", action="store_true", help="Emit the methods-summary result as JSON."
    )
    _add_manifest_argument(tree_set_methods_summary)

    tree_set_package = tree_set_subparsers.add_parser(
        "package",
        help="Build a publication-oriented uncertainty figure package for one tree set.",
    )
    tree_set_package.add_argument("tree_set", type=Path)
    tree_set_package.add_argument("--out-dir", required=True, type=Path)
    tree_set_package.add_argument("--layout", default="phylogram")
    tree_set_package.add_argument("--max-tree-count", type=int)
    tree_set_package.add_argument("--max-report-table-rows", type=int)
    tree_set_package.add_argument("--memory-warning-threshold-bytes", type=int)
    tree_set_package.add_argument(
        "--json", action="store_true", help="Emit the package result as JSON."
    )
    _add_manifest_argument(tree_set_package)

    tree_set_report = tree_set_subparsers.add_parser(
        "report",
        help="Render an HTML uncertainty report for a tree set.",
    )
    tree_set_report.add_argument("tree_set", type=Path)
    tree_set_report.add_argument("--out", required=True, type=Path)
    tree_set_report.add_argument("--max-tree-count", type=int)
    tree_set_report.add_argument("--max-report-table-rows", type=int)
    tree_set_report.add_argument("--memory-warning-threshold-bytes", type=int)
    tree_set_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(tree_set_report)


def run_tree_set_presentation_command(args: Any) -> int | None:
    if args.tree_set_command == "bootstrap-summary":
        report = write_bootstrap_tree_set_artifacts(
            args.tree_set,
            out_dir=args.out_dir,
            prefix=args.prefix,
            consensus_threshold=args.consensus_threshold,
            robust_support_threshold=args.robust_support_threshold,
        )
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=list(report.output_paths.values()),
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                warnings=report.summary_report.warnings,
                metrics={
                    "tree_count": report.summary_report.tree_count,
                    "runtime_seconds": report.summary_report.processing.runtime_seconds,
                    "peak_memory_bytes": report.summary_report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.summary_report.processing.skipped_malformed_tree_count
                    ),
                    "rooted_topology_count": (
                        report.summary_report.diversity.rooted_topology_count
                    ),
                    "unstable_branch_count": report.summary_report.unstable_branch_count,
                    "warning_count": len(report.summary_report.warnings),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "methods-summary":
        report = build_tree_set_uncertainty_method_report(args.tree_set)
        result = write_tree_set_uncertainty_methods_summary_text(args.out, report)
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=[result.output_path],
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                warnings=list(result.warnings),
                metrics={
                    "warning_count": result.warning_count,
                    "tree_count": result.report.summary.tree_count,
                    "rooted_topology_count": (
                        result.report.summary.rooted_topology_count
                    ),
                    "topology_cluster_count": result.topology_cluster_count,
                    "unstable_taxon_count": result.unstable_taxon_count,
                    "multimodal": result.report.multimodality.multimodal,
                },
                data=result,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "package":
        report = build_tree_set_uncertainty_figure_package(
            args.tree_set,
            out_dir=args.out_dir,
            layout=args.layout,
            max_tree_count=args.max_tree_count,
            max_report_table_rows=args.max_report_table_rows,
            memory_warning_threshold_bytes=args.memory_warning_threshold_bytes,
        )
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=[
                report.consensus_tree_path,
                report.consensus_figure_path,
                report.clade_support_plot_path,
                report.unstable_taxa_plot_path,
                report.topology_clusters_plot_path,
                report.unstable_taxa_table_path,
                report.topology_clusters_table_path,
                report.uncertainty_conclusions_table_path,
                report.conclusion_summary_path,
                report.legend_path,
                report.caption_path,
                report.methods_summary_path,
                report.review_path,
                report.manifest_path,
                report.reproducibility_manifest_path,
            ],
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                warnings=report.budget_report.warning_messages,
                metrics={
                    "artifact_count": 15,
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "budget_warning_count": len(report.budget_report.warning_messages),
                    "publication_ready": report.audit.publication_ready,
                    "support_labels_validated": report.audit.support_labels_validated,
                    "plotted_unstable_taxon_count": (
                        report.audit.plotted_unstable_taxon_count
                    ),
                    "plotted_topology_cluster_count": (
                        report.audit.plotted_topology_cluster_count
                    ),
                    "methods_summary_warning_count": (
                        report.methods_summary.warning_count
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command != "report":
        return None

    report = render_tree_uncertainty_report(
        tree_set_path=args.tree_set,
        out_path=args.out,
        max_tree_count=args.max_tree_count,
        max_report_table_rows=args.max_report_table_rows,
        memory_warning_threshold_bytes=args.memory_warning_threshold_bytes,
    )
    outputs = _finalize_outputs(
        args,
        command="tree-set",
        inputs=[args.tree_set],
        outputs=[args.out, report.artifact_manifest_path],
    )
    _print_result(
        build_command_result(
            command="tree-set",
            inputs=[args.tree_set],
            outputs=outputs,
            warnings=report.budget_report.warning_messages,
            metrics={
                "tree_count": report.tree_count,
                "runtime_seconds": report.processing.runtime_seconds,
                "peak_memory_bytes": report.processing.peak_memory_bytes,
                "budget_warning_count": len(report.budget_report.warning_messages),
                "methods_summary_warning_count": report.methods_summary_warning_count,
                "section_count": len(report.machine_manifest["sections"]),
                "linked_artifact_count": report.linked_artifact_count,
                "html_size_bytes": report.html_size_bytes,
                "linked_artifact_bytes": report.linked_artifact_bytes,
                "total_output_bytes": report.total_output_bytes,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
