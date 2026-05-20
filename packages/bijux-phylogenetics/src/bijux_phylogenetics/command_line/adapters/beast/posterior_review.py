from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    summarize_beast_posterior_topology_diversity,
    summarize_beast_posterior_trees,
)
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.trees.tree_sets import (
    compute_clade_frequency_table,
    compute_tree_distance_matrix,
    write_clade_frequency_table,
    write_tree_distance_matrix,
)
from bijux_phylogenetics.trees.uncertainty import (
    cluster_trees_by_topology,
    detect_unstable_clades,
    write_topology_cluster_table,
    write_unstable_clade_table,
)


def add_beast_posterior_review_commands(adapter_subparsers: Any) -> None:
    adapter_beast_consensus = adapter_subparsers.add_parser(
        "beast-consensus",
        help="Build a majority-rule consensus tree from BEAST posterior tree samples.",
    )
    adapter_beast_consensus.add_argument("input_path", type=Path)
    adapter_beast_consensus.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.25,
        help="Discard this fraction of early sampled trees before consensus building.",
    )
    adapter_beast_consensus.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Write the posterior-probability-annotated consensus tree as Newick.",
    )
    adapter_beast_consensus.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_consensus.add_argument(
        "--clade-table-out",
        type=Path,
        help="Write the retained clade-frequency ledger as TSV.",
    )
    adapter_beast_consensus.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior consensus report as JSON.",
    )
    _add_manifest_argument(adapter_beast_consensus)

    adapter_beast_diversity = adapter_subparsers.add_parser(
        "beast-diversity",
        help="Summarize topology diversity across BEAST posterior tree samples.",
    )
    adapter_beast_diversity.add_argument("input_path", type=Path)
    adapter_beast_diversity.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.25,
        help="Discard this fraction of early sampled trees before topology review.",
    )
    adapter_beast_diversity.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_diversity.add_argument(
        "--distance-out",
        type=Path,
        help="Write the pairwise RF distance ledger as TSV.",
    )
    adapter_beast_diversity.add_argument(
        "--topology-out",
        type=Path,
        help="Write the rooted topology cluster ledger as TSV.",
    )
    adapter_beast_diversity.add_argument(
        "--unstable-clade-out",
        type=Path,
        help="Write the unstable-clade ledger as TSV.",
    )
    adapter_beast_diversity.add_argument(
        "--json",
        action="store_true",
        help="Emit the topology diversity report as JSON.",
    )
    _add_manifest_argument(adapter_beast_diversity)


def run_beast_posterior_review_command(args: Any) -> int | None:
    if args.adapter_command == "beast-consensus":
        consensus_tree, report = summarize_beast_posterior_trees(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
        )
        outputs: list[Path | str] = [write_newick(args.out, consensus_tree)]
        if args.tree_set_out is not None:
            args.tree_set_out.parent.mkdir(parents=True, exist_ok=True)
            args.tree_set_out.write_text(
                report.retained_tree_set_path.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            outputs.append(args.tree_set_out)
        if args.clade_table_out is not None:
            outputs.append(
                write_clade_frequency_table(
                    args.clade_table_out,
                    compute_clade_frequency_table(report.retained_tree_set_path),
                )
            )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "total_tree_count": report.total_tree_count,
                    "kept_tree_count": report.kept_tree_count,
                    "annotated_node_count": report.annotated_node_count,
                    "clade_frequency_count": report.clade_frequency_count,
                    "burnin_fraction": report.burnin_fraction,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command != "beast-diversity":
        return None

    report = summarize_beast_posterior_topology_diversity(
        args.input_path,
        burnin_fraction=args.burnin_fraction,
    )
    outputs: list[Path | str] = []
    if args.tree_set_out is not None:
        args.tree_set_out.parent.mkdir(parents=True, exist_ok=True)
        args.tree_set_out.write_text(
            Path(report.retained_tree_set_path).read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        outputs.append(args.tree_set_out)
    if args.distance_out is not None:
        outputs.append(
            write_tree_distance_matrix(
                args.distance_out,
                compute_tree_distance_matrix(report.retained_tree_set_path),
            )
        )
    if args.topology_out is not None:
        outputs.append(
            write_topology_cluster_table(
                args.topology_out,
                cluster_trees_by_topology(report.retained_tree_set_path),
            )
        )
    if args.unstable_clade_out is not None:
        outputs.append(
            write_unstable_clade_table(
                args.unstable_clade_out,
                detect_unstable_clades(report.retained_tree_set_path),
            )
        )
    outputs = _finalize_outputs(
        args,
        command="adapter",
        inputs=[args.input_path],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="adapter",
            inputs=[args.input_path],
            outputs=outputs,
            metrics={
                "total_tree_count": report.total_tree_count,
                "kept_tree_count": report.kept_tree_count,
                "rooted_topology_count": report.rooted_topology_count,
                "dominant_topology_frequency": report.dominant_topology_frequency,
                "pair_count": report.pair_count,
                "unstable_clade_count": report.unstable_clade_count,
                "burnin_fraction": report.burnin_fraction,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
