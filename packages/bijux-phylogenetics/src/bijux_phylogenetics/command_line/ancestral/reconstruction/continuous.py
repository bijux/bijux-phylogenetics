from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.ancestral.continuous import (
    continuous_ancestral_exclusions,
    reconstruct_continuous_ancestral_states,
    summarize_continuous_ancestral_report,
    write_continuous_ancestral_exclusion_table,
    write_continuous_ancestral_summary_table,
    write_continuous_ancestral_uncertainty_table,
)
from bijux_phylogenetics.ancestral.presentation.report_rendering import (
    write_ancestral_state_table,
)
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_continuous_reconstruction_commands(ancestral_subparsers: Any) -> None:
    ancestral_continuous = ancestral_subparsers.add_parser(
        "continuous",
        help="Reconstruct ancestral states for a continuous trait.",
    )
    ancestral_continuous.add_argument("tree", type=Path)
    ancestral_continuous.add_argument("table", type=Path)
    ancestral_continuous.add_argument("--trait", required=True)
    ancestral_continuous.add_argument("--taxon-column")
    ancestral_continuous.add_argument(
        "--model", choices=("brownian", "ou"), default="brownian"
    )
    ancestral_continuous.add_argument(
        "--estimator",
        choices=("ace-pic", "anc-ml", "fast-anc", "generalized-least-squares"),
        help="Override the continuous ancestral estimator; default follows the selected model.",
    )
    ancestral_continuous.add_argument("--alpha", type=float, default=1.0)
    ancestral_continuous.add_argument("--table-out", type=Path)
    ancestral_continuous.add_argument("--summary-out", type=Path)
    ancestral_continuous.add_argument("--uncertainty-out", type=Path)
    ancestral_continuous.add_argument("--exclusions-out", type=Path)
    ancestral_continuous.add_argument(
        "--json", action="store_true", help="Emit the reconstruction as JSON."
    )
    _add_manifest_argument(ancestral_continuous)


def run_continuous_reconstruction_command(args: Any, *, parser: Any) -> int | None:
    if args.ancestral_command != "continuous":
        return None

    if args.model == "brownian" and args.estimator == "generalized-least-squares":
        parser.error(
            "continuous ancestral estimator generalized-least-squares requires model ou"
        )
    if args.model == "ou" and args.estimator in {"ace-pic", "anc-ml", "fast-anc"}:
        parser.error(
            "continuous ancestral estimators ace-pic, anc-ml, and fast-anc require model brownian"
        )
    report = reconstruct_continuous_ancestral_states(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
        estimator=args.estimator,
        alpha=args.alpha,
    )
    summary = summarize_continuous_ancestral_report(report)
    exclusions = continuous_ancestral_exclusions(report)
    outputs: list[Path | str] = []
    if args.table_out is not None:
        outputs.append(write_ancestral_state_table(args.table_out, report))
    if args.summary_out is not None:
        outputs.append(
            write_continuous_ancestral_summary_table(
                args.summary_out,
                report,
            )
        )
    if args.uncertainty_out is not None:
        outputs.append(
            write_continuous_ancestral_uncertainty_table(
                args.uncertainty_out,
                report,
            )
        )
    if args.exclusions_out is not None:
        outputs.append(
            write_continuous_ancestral_exclusion_table(
                args.exclusions_out,
                report,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="ancestral",
        inputs=[args.tree, args.table],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="ancestral",
            inputs=[args.tree, args.table],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "taxon_count": report.taxon_count,
                "estimate_count": len(report.estimates),
                "internal_node_count": summary.internal_node_count,
                "excluded_taxon_count": len(exclusions),
                "unstable_node_count": summary.unstable_node_count,
                "model": report.model,
                "estimator": report.estimator,
                "tree_is_ultrametric": summary.tree_is_ultrametric,
                "covariance_near_singular": summary.covariance_near_singular,
                "covariance_condition_number": summary.covariance_condition_number,
                "log_likelihood": summary.log_likelihood,
                "residual_sigma_squared": summary.residual_sigma_squared,
                "optimizer_name": summary.optimizer_name,
                "optimizer_converged": summary.optimizer_converged,
                "optimizer_iteration_count": summary.optimizer_iteration_count,
                "optimizer_function_evaluation_count": (
                    summary.optimizer_function_evaluation_count
                ),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
