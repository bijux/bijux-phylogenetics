from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import build_posterior_uncertainty_figure_package
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports.service import render_tree_uncertainty_report
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.trees import (
    analyze_tree_set_branch_lengths,
    cluster_trees_by_topology,
    compare_posterior_topological_diversity,
    compare_posterior_tree_sets,
    compute_clade_frequency_table,
    compute_consensus_tree,
    compute_reference_tree_clade_support,
    compute_strict_consensus_tree,
    compute_tree_distance_matrix,
    detect_posterior_topology_multimodality,
    detect_unstable_clades,
    detect_unstable_taxa,
    extract_tree_set_clades,
    load_tree_set,
    summarize_clade_credibility_conflicts,
    summarize_posterior_topology_diversity,
    summarize_tree_set_shapes,
    summarize_uncertainty_aware_conclusions,
    write_bootstrap_tree_set_artifacts,
    write_branch_length_table,
    write_clade_frequency_table,
    write_clade_table,
    write_consensus_tree,
    write_reference_tree_clade_support_table,
    write_tree_distance_distribution_table,
    write_tree_distance_matrix,
    write_tree_shape_table,
)


def add_tree_set_commands(subparsers: Any) -> None:
    tree_set = subparsers.add_parser(
        get_command_spec("tree-set").name,
        help=get_command_spec("tree-set").summary,
    )
    tree_set_subparsers = tree_set.add_subparsers(
        dest="tree_set_command", required=True
    )
    tree_set_inspect = tree_set_subparsers.add_parser(
        "inspect",
        help="Inspect a tree set for tree count and topology diversity.",
    )
    tree_set_inspect.add_argument("tree_set", type=Path)
    tree_set_inspect.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(tree_set_inspect)

    tree_set_consensus = tree_set_subparsers.add_parser(
        "consensus",
        help="Build a strict or majority-rule consensus tree from a tree set.",
    )
    tree_set_consensus.add_argument("tree_set", type=Path)
    tree_set_consensus.add_argument("--out", required=True, type=Path)
    tree_set_consensus.add_argument(
        "--method",
        choices=["majority-rule", "strict"],
        default="majority-rule",
        help="Consensus contract to apply across the governed shared taxon set.",
    )
    tree_set_consensus.add_argument(
        "--clade-frequencies-out",
        type=Path,
        help="Write the clade-frequency ledger as TSV.",
    )
    tree_set_consensus.add_argument(
        "--json", action="store_true", help="Emit the consensus report as JSON."
    )
    _add_manifest_argument(tree_set_consensus)

    tree_set_clades = tree_set_subparsers.add_parser(
        "clade-frequencies",
        help="Compute clade support frequencies across a tree set.",
    )
    tree_set_clades.add_argument("tree_set", type=Path)
    tree_set_clades.add_argument(
        "--out", type=Path, help="Write the clade-frequency table as TSV."
    )
    tree_set_clades.add_argument(
        "--json", action="store_true", help="Emit the clade-frequency report as JSON."
    )
    _add_manifest_argument(tree_set_clades)

    tree_set_support_map = tree_set_subparsers.add_parser(
        "support-map",
        help="Map tree-set clade support onto a reference tree by descendant tip set.",
    )
    tree_set_support_map.add_argument("reference_tree", type=Path)
    tree_set_support_map.add_argument("tree_set", type=Path)
    tree_set_support_map.add_argument(
        "--out", type=Path, help="Write the reference-tree support table as TSV."
    )
    tree_set_support_map.add_argument(
        "--json", action="store_true", help="Emit the support-mapping report as JSON."
    )
    _add_manifest_argument(tree_set_support_map)

    tree_set_clade_rows = tree_set_subparsers.add_parser(
        "clades",
        help="Extract one row per clade from each tree in a tree set.",
    )
    tree_set_clade_rows.add_argument("tree_set", type=Path)
    tree_set_clade_rows.add_argument("--metadata", type=Path)
    tree_set_clade_rows.add_argument("--taxon-column", type=str)
    tree_set_clade_rows.add_argument(
        "--metadata-column",
        dest="metadata_columns",
        action="append",
        default=[],
        help="Summarize this taxon-keyed metadata column for each clade.",
    )
    tree_set_clade_rows.add_argument(
        "--out", type=Path, help="Write the clade table as TSV."
    )
    tree_set_clade_rows.add_argument(
        "--json", action="store_true", help="Emit the clade table report as JSON."
    )
    _add_manifest_argument(tree_set_clade_rows)

    tree_set_shape = tree_set_subparsers.add_parser(
        "shape",
        help="Summarize tree balance, ladderization, and height across a tree set.",
    )
    tree_set_shape.add_argument("tree_set", type=Path)
    tree_set_shape.add_argument(
        "--out", type=Path, help="Write one shape-summary row per tree as TSV."
    )
    tree_set_shape.add_argument(
        "--json", action="store_true", help="Emit the tree-shape report as JSON."
    )
    _add_manifest_argument(tree_set_shape)

    tree_set_branch_lengths = tree_set_subparsers.add_parser(
        "branch-lengths",
        help="Summarize branch-length distributions across a tree set.",
    )
    tree_set_branch_lengths.add_argument("tree_set", type=Path)
    tree_set_branch_lengths.add_argument(
        "--out", type=Path, help="Write one row per branch as TSV."
    )
    tree_set_branch_lengths.add_argument(
        "--json",
        action="store_true",
        help="Emit the branch-length distribution report as JSON.",
    )
    _add_manifest_argument(tree_set_branch_lengths)

    tree_set_distances = tree_set_subparsers.add_parser(
        "distance-matrix",
        help="Compute pairwise RF distances across a tree set.",
    )
    tree_set_distances.add_argument("tree_set", type=Path)
    tree_set_distances.add_argument(
        "--out", type=Path, help="Write the pairwise distance table as TSV."
    )
    tree_set_distances.add_argument(
        "--json", action="store_true", help="Emit the distance report as JSON."
    )
    _add_manifest_argument(tree_set_distances)

    tree_set_diversity_summary = tree_set_subparsers.add_parser(
        "diversity",
        help="Summarize topology diversity and RF distribution across one tree set.",
    )
    tree_set_diversity_summary.add_argument("tree_set", type=Path)
    tree_set_diversity_summary.add_argument(
        "--out",
        type=Path,
        help="Write the RF-distribution table as TSV.",
    )
    tree_set_diversity_summary.add_argument(
        "--json", action="store_true", help="Emit the diversity report as JSON."
    )
    _add_manifest_argument(tree_set_diversity_summary)

    tree_set_clusters = tree_set_subparsers.add_parser(
        "cluster",
        help="Cluster trees by identical rooted topology signatures.",
    )
    tree_set_clusters.add_argument("tree_set", type=Path)
    tree_set_clusters.add_argument(
        "--json", action="store_true", help="Emit the cluster report as JSON."
    )
    _add_manifest_argument(tree_set_clusters)

    tree_set_unstable_taxa = tree_set_subparsers.add_parser(
        "unstable-taxa",
        help="Detect taxa with inconsistent placements across a tree set.",
    )
    tree_set_unstable_taxa.add_argument("tree_set", type=Path)
    tree_set_unstable_taxa.add_argument(
        "--json", action="store_true", help="Emit the instability report as JSON."
    )
    _add_manifest_argument(tree_set_unstable_taxa)

    tree_set_unstable_clades = tree_set_subparsers.add_parser(
        "unstable-clades",
        help="Detect non-unanimous and conflicting clades across a tree set.",
    )
    tree_set_unstable_clades.add_argument("tree_set", type=Path)
    tree_set_unstable_clades.add_argument(
        "--json", action="store_true", help="Emit the instability report as JSON."
    )
    _add_manifest_argument(tree_set_unstable_clades)

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

    tree_set_compare = tree_set_subparsers.add_parser(
        "compare",
        help="Compare two posterior tree sets over clade support and topology distance.",
    )
    tree_set_compare.add_argument("left", type=Path)
    tree_set_compare.add_argument("right", type=Path)
    tree_set_compare.add_argument(
        "--json", action="store_true", help="Emit the comparison report as JSON."
    )
    _add_manifest_argument(tree_set_compare)

    tree_set_diversity = tree_set_subparsers.add_parser(
        "diversity-compare",
        help="Compare posterior topological diversity across two analyses.",
    )
    tree_set_diversity.add_argument("left", type=Path)
    tree_set_diversity.add_argument("right", type=Path)
    tree_set_diversity.add_argument(
        "--json", action="store_true", help="Emit the diversity report as JSON."
    )
    _add_manifest_argument(tree_set_diversity)

    tree_set_multimodality = tree_set_subparsers.add_parser(
        "multimodality",
        help="Detect multimodal posterior topology distributions.",
    )
    tree_set_multimodality.add_argument("tree_set", type=Path)
    tree_set_multimodality.add_argument("--min-mode-frequency", type=float, default=0.2)
    tree_set_multimodality.add_argument("--min-mode-count", type=int, default=2)
    tree_set_multimodality.add_argument(
        "--json", action="store_true", help="Emit the multimodality report as JSON."
    )
    _add_manifest_argument(tree_set_multimodality)

    tree_set_conflicts = tree_set_subparsers.add_parser(
        "clade-conflicts",
        help="Summarize conflicting high-credibility clades across a posterior tree set.",
    )
    tree_set_conflicts.add_argument("tree_set", type=Path)
    tree_set_conflicts.add_argument("--credibility-threshold", type=float, default=0.5)
    tree_set_conflicts.add_argument(
        "--json", action="store_true", help="Emit the clade-conflict report as JSON."
    )
    _add_manifest_argument(tree_set_conflicts)

    tree_set_summary = tree_set_subparsers.add_parser(
        "conclusion-summary",
        help="Summarize robust, uncertain, and conflict-prone clades from posterior uncertainty.",
    )
    tree_set_summary.add_argument("tree_set", type=Path)
    tree_set_summary.add_argument("--robust-threshold", type=float, default=0.9)
    tree_set_summary.add_argument("--uncertain-min-frequency", type=float, default=0.3)
    tree_set_summary.add_argument("--uncertain-max-frequency", type=float, default=0.7)
    tree_set_summary.add_argument("--credibility-threshold", type=float, default=0.5)
    tree_set_summary.add_argument(
        "--json", action="store_true", help="Emit the conclusion summary as JSON."
    )
    _add_manifest_argument(tree_set_summary)

    tree_set_package = tree_set_subparsers.add_parser(
        "package",
        help="Build a posterior uncertainty figure package for one tree set.",
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


def run_tree_set_command(args: Any) -> int:
    if args.tree_set_command == "inspect":
        report = load_tree_set(args.tree_set)
        outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "shared_taxon_count": len(report.shared_taxa),
                    "rooted_topology_count": report.rooted_topology_count,
                    "unrooted_topology_count": report.unrooted_topology_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "consensus":
        if args.method == "strict":
            tree, report = compute_strict_consensus_tree(args.tree_set)
        else:
            tree, report = compute_consensus_tree(args.tree_set)
        output_paths = [write_consensus_tree(args.out, tree)]
        if args.clade_frequencies_out is not None:
            frequency_report = compute_clade_frequency_table(args.tree_set)
            output_paths.append(
                write_clade_frequency_table(
                    args.clade_frequencies_out,
                    frequency_report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "shared_taxon_count": len(report.shared_taxa),
                    "consensus_method": report.consensus_method,
                    "consensus_threshold": report.consensus_threshold,
                    "included_clade_count": report.included_clade_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "clade-frequencies":
        report = compute_clade_frequency_table(args.tree_set)
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_clade_frequency_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "clade_count": len(report.clade_frequencies),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "support-map":
        report = compute_reference_tree_clade_support(
            args.reference_tree,
            args.tree_set,
        )
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_reference_tree_clade_support_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.reference_tree, args.tree_set],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.reference_tree, args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "shared_taxon_count": len(report.shared_taxa),
                    "supported_clade_count": report.supported_clade_count,
                    "absent_clade_count": report.absent_clade_count,
                    "unscored_clade_count": report.unscored_clade_count,
                    "reference_internal_node_count": len(report.rows),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "clades":
        report = extract_tree_set_clades(
            args.tree_set,
            metadata_path=args.metadata,
            taxon_column=args.taxon_column,
            metadata_columns=args.metadata_columns or None,
        )
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_clade_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "clade_count": len(report.rows),
                    "metadata_column_count": len(report.metadata_columns),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "shape":
        report = summarize_tree_set_shapes(args.tree_set)
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_tree_shape_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "balanced_tree_count": report.aggregate.balanced_tree_count,
                    "ladderized_tree_count": report.aggregate.ladderized_tree_count,
                    "star_like_tree_count": report.aggregate.star_like_tree_count,
                    "comb_like_tree_count": report.aggregate.comb_like_tree_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "branch-lengths":
        report = analyze_tree_set_branch_lengths(args.tree_set)
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_branch_length_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "branch_count": report.aggregate.branch_count,
                    "zero_length_branch_count": report.aggregate.zero_length_branch_count,
                    "negative_branch_count": report.aggregate.negative_branch_count,
                    "long_outlier_count": report.aggregate.long_outlier_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "distance-matrix":
        report = compute_tree_distance_matrix(args.tree_set)
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_tree_distance_matrix(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "pair_count": len(report.pairs),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "diversity":
        report = summarize_posterior_topology_diversity(args.tree_set)
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_tree_distance_distribution_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "rooted_topology_count": report.rooted_topology_count,
                    "pair_count": report.pair_count,
                    "rf_bucket_count": len(report.rf_distribution),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "cluster":
        report = cluster_trees_by_topology(args.tree_set)
        outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "cluster_count": len(report.clusters),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "unstable-taxa":
        report = detect_unstable_taxa(args.tree_set)
        outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "unstable_taxon_count": len(report.taxa),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "unstable-clades":
        report = detect_unstable_clades(args.tree_set)
        outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "unstable_clade_count": len(report.clades),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

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

    if args.tree_set_command == "compare":
        report = compare_posterior_tree_sets(args.left, args.right)
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.left, args.right],
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.left, args.right],
                outputs=outputs,
                metrics={
                    "left_tree_count": report.left_tree_count,
                    "right_tree_count": report.right_tree_count,
                    "shared_rooted_topology_count": report.shared_rooted_topology_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "diversity-compare":
        report = compare_posterior_topological_diversity(args.left, args.right)
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.left, args.right],
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.left, args.right],
                outputs=outputs,
                metrics={
                    "left_rooted_topology_count": (
                        report.left_summary.rooted_topology_count
                    ),
                    "right_rooted_topology_count": (
                        report.right_summary.rooted_topology_count
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "multimodality":
        report = detect_posterior_topology_multimodality(
            args.tree_set,
            min_mode_frequency=args.min_mode_frequency,
            min_mode_count=args.min_mode_count,
        )
        outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "mode_count": report.mode_count,
                    "multimodal": report.multimodal,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "clade-conflicts":
        report = summarize_clade_credibility_conflicts(
            args.tree_set,
            credibility_threshold=args.credibility_threshold,
        )
        outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={"conflict_count": report.conflict_count},
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "conclusion-summary":
        report = summarize_uncertainty_aware_conclusions(
            args.tree_set,
            robust_threshold=args.robust_threshold,
            uncertain_min_frequency=args.uncertain_min_frequency,
            uncertain_max_frequency=args.uncertain_max_frequency,
            credibility_threshold=args.credibility_threshold,
        )
        outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "robust_clade_count": report.robust_clade_count,
                    "uncertain_clade_count": report.uncertain_clade_count,
                    "conflicting_clade_count": report.conflicting_clade_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "package":
        report = build_posterior_uncertainty_figure_package(
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
                report.clade_frequency_plot_path,
                report.unstable_taxa_table_path,
                report.topology_clusters_table_path,
                report.conclusion_summary_path,
                report.manifest_path,
            ],
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                warnings=report.budget_report.warning_messages,
                metrics={
                    "artifact_count": 7,
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "budget_warning_count": len(report.budget_report.warning_messages),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

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
