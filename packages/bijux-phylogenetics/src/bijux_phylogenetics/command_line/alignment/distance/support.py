from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_distance_tree_method_argument,
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.distance import (
    bootstrap_distance_trees,
    summarize_distance_bootstrap_support,
    write_distance_bootstrap_draws,
    write_distance_bootstrap_support,
)
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.simulation import write_tree_set

from .shared import (
    add_ambiguity_policy_option,
    add_distance_model_option,
    add_gap_handling_option,
)


def add_distance_support_commands(alignment_subparsers: Any) -> None:
    alignment_bootstrap_tree = alignment_subparsers.add_parser(
        "bootstrap-tree",
        help="Bootstrap a distance tree by resampling alignment sites with replacement.",
    )
    alignment_bootstrap_tree.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_bootstrap_tree)
    add_distance_model_option(alignment_bootstrap_tree)
    add_gap_handling_option(alignment_bootstrap_tree)
    add_ambiguity_policy_option(alignment_bootstrap_tree)
    alignment_bootstrap_tree.add_argument("--replicates", type=int, default=100)
    alignment_bootstrap_tree.add_argument("--seed", type=int, default=1)
    alignment_bootstrap_tree.add_argument(
        "--support-out", type=Path, help="Write bootstrap clade support as TSV."
    )
    alignment_bootstrap_tree.add_argument(
        "--tree-set-out", type=Path, help="Write bootstrap replicate trees as Newick."
    )
    alignment_bootstrap_tree.add_argument(
        "--draws-out",
        type=Path,
        help="Write the deterministic bootstrap site-draw ledger as TSV.",
    )
    alignment_bootstrap_tree.add_argument(
        "--json", action="store_true", help="Emit the bootstrap report as JSON."
    )
    _add_manifest_argument(alignment_bootstrap_tree)

    alignment_bootstrap_summary = alignment_subparsers.add_parser(
        "distance-support-summary",
        help="Summarize consensus clade support across distance-bootstrap replicates.",
    )
    alignment_bootstrap_summary.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_bootstrap_summary)
    add_distance_model_option(alignment_bootstrap_summary)
    add_gap_handling_option(alignment_bootstrap_summary)
    add_ambiguity_policy_option(alignment_bootstrap_summary)
    alignment_bootstrap_summary.add_argument("--replicates", type=int, default=25)
    alignment_bootstrap_summary.add_argument("--seed", type=int, default=1)
    alignment_bootstrap_summary.add_argument(
        "--json", action="store_true", help="Emit the support summary as JSON."
    )
    _add_manifest_argument(alignment_bootstrap_summary)


def run_distance_support_command(args: Any) -> int | None:
    if args.alignment_command == "bootstrap-tree":
        trees, report = bootstrap_distance_trees(
            args.alignment,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
            replicates=args.replicates,
            seed=args.seed,
        )
        outputs: list[Path | str] = []
        if args.support_out is not None:
            outputs.append(write_distance_bootstrap_support(args.support_out, report))
        if args.tree_set_out is not None:
            outputs.append(write_tree_set(args.tree_set_out, trees))
        if args.draws_out is not None:
            outputs.append(write_distance_bootstrap_draws(args.draws_out, report))
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "replicate_count": report.tree_count,
                    "support_row_count": len(report.support),
                    "method": report.method,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.alignment_command == "distance-support-summary":
        report = summarize_distance_bootstrap_support(
            bootstrap_distance_trees(
                args.alignment,
                method=args.method,
                model=args.model,
                gap_handling=args.gap_handling,
                ambiguity_policy=args.ambiguity_policy,
                replicates=args.replicates,
                seed=args.seed,
            )[1]
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "clade_count": report.clade_count,
                    "weak_clade_count": report.weak_clade_count,
                    "replicates": report.replicates,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
