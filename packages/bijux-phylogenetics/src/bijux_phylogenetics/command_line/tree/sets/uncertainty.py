from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.trees import (
    detect_posterior_topology_multimodality,
    detect_rogue_taxa,
    detect_unstable_clades,
    detect_unstable_taxa,
    summarize_clade_credibility_conflicts,
    summarize_uncertainty_aware_conclusions,
    write_consensus_tree,
    write_gene_tree_conflict_artifacts,
    write_rogue_taxon_table,
)


def add_tree_set_uncertainty_commands(tree_set_subparsers: Any) -> None:
    tree_set_rogue_taxa = tree_set_subparsers.add_parser(
        "rogue-taxa",
        help="Rank taxa by the consensus-resolution, support, and RF-stability improvement gained by removing them.",
    )
    tree_set_rogue_taxa.add_argument("tree_set", type=Path)
    tree_set_rogue_taxa.add_argument("--out-dir", required=True, type=Path)
    tree_set_rogue_taxa.add_argument(
        "--consensus-threshold",
        type=float,
        default=0.5,
        help="Clade-frequency threshold used when scoring baseline and pruned consensus trees.",
    )
    tree_set_rogue_taxa.add_argument(
        "--json", action="store_true", help="Emit the rogue-taxon report as JSON."
    )
    _add_manifest_argument(tree_set_rogue_taxa)

    tree_set_gene_tree_conflicts = tree_set_subparsers.add_parser(
        "gene-tree-conflicts",
        help="Write clade, quartet, rogue-taxon, and conflict tables for one gene-tree set.",
    )
    tree_set_gene_tree_conflicts.add_argument("tree_set", type=Path)
    tree_set_gene_tree_conflicts.add_argument("--out-dir", required=True, type=Path)
    tree_set_gene_tree_conflicts.add_argument(
        "--prefix",
        default="gene-tree-conflicts",
        help="Prefix for written artifacts.",
    )
    tree_set_gene_tree_conflicts.add_argument(
        "--credibility-threshold",
        type=float,
        default=0.5,
        help="Clade-frequency threshold used when flagging conflicting high-credibility clades.",
    )
    tree_set_gene_tree_conflicts.add_argument(
        "--rogue-consensus-threshold",
        type=float,
        default=0.5,
        help="Consensus threshold used when ranking rogue taxa.",
    )
    tree_set_gene_tree_conflicts.add_argument(
        "--json",
        action="store_true",
        help="Emit the gene-tree conflict bundle report as JSON.",
    )
    _add_manifest_argument(tree_set_gene_tree_conflicts)

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


def run_tree_set_uncertainty_command(args: Any) -> int | None:
    if args.tree_set_command == "gene-tree-conflicts":
        report = write_gene_tree_conflict_artifacts(
            args.tree_set,
            out_dir=args.out_dir,
            prefix=args.prefix,
            credibility_threshold=args.credibility_threshold,
            rogue_consensus_threshold=args.rogue_consensus_threshold,
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
                metrics={
                    "tree_count": report.summary_report.tree_count,
                    "runtime_seconds": report.summary_report.processing.runtime_seconds,
                    "peak_memory_bytes": report.summary_report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.summary_report.processing.skipped_malformed_tree_count
                    ),
                    "shared_taxon_count": len(report.summary_report.shared_taxa),
                    "reference_tree_frequency": (
                        report.summary_report.reference_tree.frequency
                    ),
                    "clade_count": len(
                        report.summary_report.clade_frequencies.clade_frequencies
                    ),
                    "quartet_branch_count": (
                        report.summary_report.quartet_concordance.branch_count
                    ),
                    "conflict_count": report.summary_report.clade_conflicts.conflict_count,
                    "rogue_taxon_count": len(report.summary_report.rogue_taxa.rows),
                    "top_ranked_rogue_taxon": (
                        report.summary_report.rogue_taxa.rows[0].taxon
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "rogue-taxa":
        report = detect_rogue_taxa(
            args.tree_set,
            consensus_threshold=args.consensus_threshold,
        )
        args.out_dir.mkdir(parents=True, exist_ok=True)
        ranking_out = write_rogue_taxon_table(
            args.out_dir / "rogue-taxon-ranking.tsv",
            report,
        )
        baseline_consensus_out = write_consensus_tree(
            args.out_dir / "baseline-consensus.nwk",
            loads_newick(report.baseline_consensus_newick),
        )
        best_consensus_out = write_consensus_tree(
            args.out_dir / "best-rogue-taxon-consensus.nwk",
            loads_newick(report.rows[0].pruned_consensus_newick),
        )
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=[ranking_out, baseline_consensus_out, best_consensus_out],
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
                    "candidate_taxon_count": len(report.rows),
                    "top_ranked_taxon": report.rows[0].taxon,
                    "consensus_threshold": report.consensus_threshold,
                    "baseline_consensus_resolution": (
                        report.baseline_consensus_resolution
                    ),
                    "top_ranked_consensus_resolution": (
                        report.rows[0].pruned_consensus_resolution
                    ),
                    "baseline_mean_support_percent": (
                        report.baseline_mean_support_percent
                    ),
                    "top_ranked_mean_support_percent": (
                        report.rows[0].pruned_mean_support_percent
                    ),
                    "baseline_mean_normalized_robinson_foulds": (
                        report.baseline_mean_normalized_robinson_foulds
                    ),
                    "top_ranked_mean_normalized_robinson_foulds": (
                        report.rows[0].pruned_mean_normalized_robinson_foulds
                    ),
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

    return None
