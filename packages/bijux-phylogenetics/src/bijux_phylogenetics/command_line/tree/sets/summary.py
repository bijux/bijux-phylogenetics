from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.trees import (
    build_quartet_puzzling_consensus,
    compute_candidate_tree_quartet_score,
    compute_clade_compatibility_graph,
    compute_clade_frequency_table,
    compute_consensus_tree,
    compute_credible_clade_set,
    compute_gene_tree_quartet_concordance_factors,
    compute_majority_rule_extended_consensus,
    compute_maximum_clade_credibility_tree,
    compute_posterior_clade_correlation_matrix,
    compute_posterior_tree_distance_diagnostics,
    compute_reference_tree_clade_support,
    compute_reference_tree_quartet_support,
    compute_strict_consensus_tree,
    extract_tree_set_clades,
    load_tree_set,
    summarize_posterior_agreement_subtree,
    summarize_posterior_branch_lengths,
    summarize_posterior_node_ages,
    write_candidate_tree_quartet_score_table,
    write_clade_compatibility_edge_table,
    write_clade_compatibility_graph_dot,
    write_clade_compatibility_node_table,
    write_clade_frequency_table,
    write_clade_table,
    write_consensus_tree,
    write_credible_clade_set_artifacts,
    write_gene_tree_quartet_concordance_table,
    write_majority_rule_extended_consensus_artifacts,
    write_maximum_clade_credibility_artifacts,
    write_posterior_agreement_subtree_artifacts,
    write_posterior_branch_length_summary_table,
    write_posterior_clade_correlation_artifacts,
    write_posterior_node_age_summary_table,
    write_posterior_tree_distance_artifacts,
    write_quartet_puzzling_artifacts,
    write_reference_tree_clade_support_table,
    write_reference_tree_quartet_support_table,
)


def add_tree_set_summary_commands(tree_set_subparsers: Any) -> None:
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

    tree_set_compatibility_graph = tree_set_subparsers.add_parser(
        "compatibility-graph",
        help="Build a clade compatibility graph where nodes are informative clades and edges are compatibility or conflict relationships.",
    )
    tree_set_compatibility_graph.add_argument("tree_set", type=Path)
    tree_set_compatibility_graph.add_argument("--out-dir", required=True, type=Path)
    tree_set_compatibility_graph.add_argument(
        "--json",
        action="store_true",
        help="Emit the compatibility graph report as JSON.",
    )
    _add_manifest_argument(tree_set_compatibility_graph)

    tree_set_extended_consensus = tree_set_subparsers.add_parser(
        "majority-rule-extended-consensus",
        help="Build a majority-rule extended consensus tree and emit the accepted inclusion order plus rejected conflict ledger.",
    )
    tree_set_extended_consensus.add_argument("tree_set", type=Path)
    tree_set_extended_consensus.add_argument("--out-dir", required=True, type=Path)
    tree_set_extended_consensus.add_argument(
        "--json",
        action="store_true",
        help="Emit the extended consensus report as JSON.",
    )
    _add_manifest_argument(tree_set_extended_consensus)

    tree_set_maximum_clade_credibility = tree_set_subparsers.add_parser(
        "maximum-clade-credibility",
        help="Select the sampled tree with the highest summed posterior clade credibility and emit the candidate score ledger.",
    )
    tree_set_maximum_clade_credibility.add_argument("tree_set", type=Path)
    tree_set_maximum_clade_credibility.add_argument(
        "--out-dir", required=True, type=Path
    )
    tree_set_maximum_clade_credibility.add_argument(
        "--json",
        action="store_true",
        help="Emit the maximum-clade-credibility report as JSON.",
    )
    _add_manifest_argument(tree_set_maximum_clade_credibility)

    tree_set_credible_clade_set = tree_set_subparsers.add_parser(
        "credible-clade-set",
        help="Select clades by descending posterior frequency until a credible-threshold cutoff is reached, and emit included plus excluded ledgers.",
    )
    tree_set_credible_clade_set.add_argument("tree_set", type=Path)
    tree_set_credible_clade_set.add_argument("--out-dir", required=True, type=Path)
    tree_set_credible_clade_set.add_argument(
        "--credible-threshold",
        type=float,
        default=0.95,
        help="Posterior clade-frequency mass that the included clade set must reach or exceed.",
    )
    tree_set_credible_clade_set.add_argument(
        "--json",
        action="store_true",
        help="Emit the credible-clade-set report as JSON.",
    )
    _add_manifest_argument(tree_set_credible_clade_set)

    tree_set_posterior_branch_lengths = tree_set_subparsers.add_parser(
        "posterior-branch-lengths",
        help="Summarize posterior branch lengths by rooted clade identity across one tree set.",
    )
    tree_set_posterior_branch_lengths.add_argument("tree_set", type=Path)
    tree_set_posterior_branch_lengths.add_argument(
        "--out",
        type=Path,
        help="Write the posterior branch-length summary table as TSV.",
    )
    tree_set_posterior_branch_lengths.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior branch-length summary report as JSON.",
    )
    _add_manifest_argument(tree_set_posterior_branch_lengths)

    tree_set_posterior_node_ages = tree_set_subparsers.add_parser(
        "posterior-node-ages",
        help="Summarize posterior node ages by rooted clade identity across one dated tree set.",
    )
    tree_set_posterior_node_ages.add_argument("tree_set", type=Path)
    tree_set_posterior_node_ages.add_argument(
        "--out",
        type=Path,
        help="Write the posterior node-age summary table as TSV.",
    )
    tree_set_posterior_node_ages.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior node-age summary report as JSON.",
    )
    _add_manifest_argument(tree_set_posterior_node_ages)

    tree_set_posterior_agreement_subtree = tree_set_subparsers.add_parser(
        "posterior-agreement-subtree",
        help="Find the largest retained taxon subset whose pruned posterior trees share one rooted topology.",
    )
    tree_set_posterior_agreement_subtree.add_argument("tree_set", type=Path)
    tree_set_posterior_agreement_subtree.add_argument(
        "--out-dir", required=True, type=Path
    )
    tree_set_posterior_agreement_subtree.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior agreement-subtree report as JSON.",
    )
    _add_manifest_argument(tree_set_posterior_agreement_subtree)

    tree_set_posterior_clade_correlation = tree_set_subparsers.add_parser(
        "posterior-clade-correlation-matrix",
        help="Compute informative clade co-occurrence and binary correlation across one posterior tree set.",
    )
    tree_set_posterior_clade_correlation.add_argument("tree_set", type=Path)
    tree_set_posterior_clade_correlation.add_argument(
        "--out-dir", required=True, type=Path
    )
    tree_set_posterior_clade_correlation.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior clade-correlation report as JSON.",
    )
    _add_manifest_argument(tree_set_posterior_clade_correlation)

    tree_set_posterior_tree_distances = tree_set_subparsers.add_parser(
        "posterior-tree-distance-diagnostics",
        help="Compare every posterior tree against MCC and consensus references by RF and branch-score distance.",
    )
    tree_set_posterior_tree_distances.add_argument("tree_set", type=Path)
    tree_set_posterior_tree_distances.add_argument(
        "--out-dir", required=True, type=Path
    )
    tree_set_posterior_tree_distances.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior tree-distance diagnostics report as JSON.",
    )
    _add_manifest_argument(tree_set_posterior_tree_distances)

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

    tree_set_quartet_support = tree_set_subparsers.add_parser(
        "quartet-support",
        help="Score each informative reference-tree branch by concordant, discordant, and uninformative induced quartets across a tree set.",
    )
    tree_set_quartet_support.add_argument("reference_tree", type=Path)
    tree_set_quartet_support.add_argument("tree_set", type=Path)
    tree_set_quartet_support.add_argument(
        "--out",
        type=Path,
        help="Write the reference-tree quartet support table as TSV.",
    )
    tree_set_quartet_support.add_argument(
        "--json", action="store_true", help="Emit the quartet-support report as JSON."
    )
    _add_manifest_argument(tree_set_quartet_support)

    tree_set_quartet_concordance = tree_set_subparsers.add_parser(
        "quartet-concordance-factors",
        help="Score each informative species-tree branch by concordant, discordant-first, discordant-second, and uninformative quartet votes across a gene-tree set.",
    )
    tree_set_quartet_concordance.add_argument("reference_tree", type=Path)
    tree_set_quartet_concordance.add_argument("tree_set", type=Path)
    tree_set_quartet_concordance.add_argument(
        "--out",
        type=Path,
        help="Write the species-tree quartet concordance-factor table as TSV.",
    )
    tree_set_quartet_concordance.add_argument(
        "--json",
        action="store_true",
        help="Emit the quartet-concordance report as JSON.",
    )
    _add_manifest_argument(tree_set_quartet_concordance)

    tree_set_quartet_score = tree_set_subparsers.add_parser(
        "quartet-score",
        help="Score one candidate species tree by quartet agreement across a gene-tree set.",
    )
    tree_set_quartet_score.add_argument("reference_tree", type=Path)
    tree_set_quartet_score.add_argument("tree_set", type=Path)
    tree_set_quartet_score.add_argument(
        "--out",
        type=Path,
        help="Write the candidate-tree quartet score table as TSV.",
    )
    tree_set_quartet_score.add_argument(
        "--json",
        action="store_true",
        help="Emit the quartet-score report as JSON.",
    )
    _add_manifest_argument(tree_set_quartet_score)

    tree_set_quartet_puzzling = tree_set_subparsers.add_parser(
        "quartet-puzzling",
        help="Assemble quartet-scored trees across deterministic taxon orders and summarize a consensus tree.",
    )
    tree_set_quartet_puzzling.add_argument("tree_set", type=Path)
    tree_set_quartet_puzzling.add_argument("--out-dir", required=True, type=Path)
    tree_set_quartet_puzzling.add_argument(
        "--max-order-count",
        type=int,
        default=32,
        help="Cap the number of deterministic taxon orders assembled into the consensus tree.",
    )
    tree_set_quartet_puzzling.add_argument(
        "--random-seed",
        type=int,
        default=0,
        help="Seed used when additional shuffled taxon orders are needed beyond the deterministic base set.",
    )
    tree_set_quartet_puzzling.add_argument(
        "--consensus-threshold",
        type=float,
        default=0.5,
        help="Clade-frequency threshold applied when summarizing the assembled tree set.",
    )
    tree_set_quartet_puzzling.add_argument(
        "--json", action="store_true", help="Emit the quartet-puzzling report as JSON."
    )
    _add_manifest_argument(tree_set_quartet_puzzling)

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


def run_tree_set_summary_command(args: Any) -> int | None:
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

    if args.tree_set_command == "compatibility-graph":
        report = compute_clade_compatibility_graph(args.tree_set)
        args.out_dir.mkdir(parents=True, exist_ok=True)
        node_out = write_clade_compatibility_node_table(
            args.out_dir / "clade-compatibility-nodes.tsv",
            report,
        )
        edge_out = write_clade_compatibility_edge_table(
            args.out_dir / "clade-compatibility-edges.tsv",
            report,
        )
        dot_out = write_clade_compatibility_graph_dot(
            args.out_dir / "clade-compatibility.dot",
            report,
        )
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=[node_out, edge_out, dot_out],
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
                    "node_count": report.node_count,
                    "edge_count": report.edge_count,
                    "compatible_edge_count": report.compatible_edge_count,
                    "conflict_edge_count": report.conflict_edge_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "majority-rule-extended-consensus":
        report = compute_majority_rule_extended_consensus(args.tree_set)[1]
        output_paths = write_majority_rule_extended_consensus_artifacts(
            args.out_dir,
            report,
        )
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=list(output_paths.values()),
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
                    "included_clade_count": report.included_clade_count,
                    "majority_included_clade_count": (
                        report.majority_included_clade_count
                    ),
                    "extension_included_clade_count": (
                        report.extension_included_clade_count
                    ),
                    "rejected_conflict_count": report.rejected_conflict_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "maximum-clade-credibility":
        tree, report = compute_maximum_clade_credibility_tree(args.tree_set)
        output_paths = write_maximum_clade_credibility_artifacts(
            args.out_dir,
            tree,
            report,
        )
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=list(output_paths.values()),
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
                    "rooted_topology_count": report.rooted_topology_count,
                    "selected_tree_index": report.selected_tree_index,
                    "clade_credibility_score": report.clade_credibility_score,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "credible-clade-set":
        report = compute_credible_clade_set(
            args.tree_set,
            credible_threshold=args.credible_threshold,
        )
        output_paths = write_credible_clade_set_artifacts(args.out_dir, report)
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=list(output_paths.values()),
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
                    "credible_threshold": report.credible_threshold,
                    "included_clade_count": report.included_clade_count,
                    "excluded_clade_count": report.excluded_clade_count,
                    "included_cumulative_frequency": (
                        report.included_cumulative_frequency
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "posterior-branch-lengths":
        report = summarize_posterior_branch_lengths(args.tree_set)
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(
                write_posterior_branch_length_summary_table(args.out, report)
            )
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
                    "shared_taxon_count": len(report.shared_taxa),
                    "hpd_mass": report.hpd_mass,
                    "branch_summary_count": len(report.rows),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "posterior-node-ages":
        report = summarize_posterior_node_ages(args.tree_set)
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_posterior_node_age_summary_table(args.out, report))
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
                    "shared_taxon_count": len(report.shared_taxa),
                    "hpd_mass": report.hpd_mass,
                    "node_age_summary_count": len(report.rows),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "posterior-agreement-subtree":
        tree, report = summarize_posterior_agreement_subtree(args.tree_set)
        output_paths = write_posterior_agreement_subtree_artifacts(
            args.out_dir,
            tree,
            report,
        )
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=list(output_paths.values()),
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
                    "evaluated_candidate_count": report.evaluated_candidate_count,
                    "retained_taxon_count": len(report.retained_taxa),
                    "agreement_removed_taxon_count": len(report.agreement_removed_taxa),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "posterior-clade-correlation-matrix":
        report = compute_posterior_clade_correlation_matrix(args.tree_set)
        output_paths = write_posterior_clade_correlation_artifacts(
            args.out_dir,
            report,
        )
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=list(output_paths.values()),
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
                    "clade_count": report.clade_count,
                    "pair_count": report.pair_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "posterior-tree-distance-diagnostics":
        report = compute_posterior_tree_distance_diagnostics(args.tree_set)
        mcc_tree, _mcc_report = compute_maximum_clade_credibility_tree(args.tree_set)
        consensus_tree, _consensus_report = compute_consensus_tree(args.tree_set)
        output_paths = write_posterior_tree_distance_artifacts(
            args.out_dir,
            mcc_tree,
            consensus_tree,
            report,
        )
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=list(output_paths.values()),
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
                    "row_count": report.row_count,
                    "distribution_row_count": report.distribution_row_count,
                    "maximum_clade_credibility_tree_index": (
                        report.maximum_clade_credibility_tree_index
                    ),
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

    if args.tree_set_command == "quartet-support":
        report = compute_reference_tree_quartet_support(
            args.reference_tree,
            args.tree_set,
        )
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_reference_tree_quartet_support_table(args.out, report))
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
                    "branch_count": report.branch_count,
                    "total_quartet_count": report.total_quartet_count,
                    "concordant_quartet_count": report.concordant_quartet_count,
                    "discordant_quartet_count": report.discordant_quartet_count,
                    "uninformative_quartet_count": report.uninformative_quartet_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "quartet-concordance-factors":
        report = compute_gene_tree_quartet_concordance_factors(
            args.reference_tree,
            args.tree_set,
        )
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_gene_tree_quartet_concordance_table(args.out, report))
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
                    "branch_count": report.branch_count,
                    "total_quartet_count": report.total_quartet_count,
                    "concordant_quartet_count": report.concordant_quartet_count,
                    "discordant_first_quartet_count": (
                        report.discordant_first_quartet_count
                    ),
                    "discordant_second_quartet_count": (
                        report.discordant_second_quartet_count
                    ),
                    "uninformative_quartet_count": report.uninformative_quartet_count,
                    "informative_quartet_count": report.informative_quartet_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "quartet-score":
        report = compute_candidate_tree_quartet_score(
            args.reference_tree,
            args.tree_set,
        )
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_candidate_tree_quartet_score_table(args.out, report))
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
                    "branch_count": report.branch_count,
                    "total_quartet_count": report.total_quartet_count,
                    "concordant_quartet_count": report.concordant_quartet_count,
                    "discordant_first_quartet_count": (
                        report.discordant_first_quartet_count
                    ),
                    "discordant_second_quartet_count": (
                        report.discordant_second_quartet_count
                    ),
                    "uninformative_quartet_count": report.uninformative_quartet_count,
                    "informative_quartet_count": report.informative_quartet_count,
                    "quartet_score": report.quartet_score,
                    "normalized_quartet_score": report.normalized_quartet_score,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "quartet-puzzling":
        _tree, report = build_quartet_puzzling_consensus(
            args.tree_set,
            max_order_count=args.max_order_count,
            random_seed=args.random_seed,
            consensus_threshold=args.consensus_threshold,
        )
        artifact_paths = write_quartet_puzzling_artifacts(args.out_dir, report)
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=list(artifact_paths.values()),
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
                    "quartet_count": report.quartet_count,
                    "assembly_count": report.assembly_count,
                    "unique_assembled_topology_count": (
                        report.unique_assembled_topology_count
                    ),
                    "canonical_root_taxon": report.canonical_root_taxon,
                    "consensus_threshold": report.consensus_threshold,
                    "included_clade_count": report.included_clade_count,
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

    return None
