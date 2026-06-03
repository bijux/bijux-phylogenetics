from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.signal import (
    summarize_phylogenetic_signal,
    write_phylogenetic_signal_permutation_table,
    write_phylogenetic_signal_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_phylogenetic_signal_commands(comparative_subparsers: Any) -> None:
    comparative_signal = comparative_subparsers.add_parser(
        "signal",
        help="Estimate phylogenetic signal metrics for one numeric trait.",
    )
    comparative_signal.add_argument("tree", type=Path)
    comparative_signal.add_argument("table", type=Path)
    comparative_signal.add_argument("--trait", required=True)
    comparative_signal.add_argument("--taxon-column")
    comparative_signal.add_argument(
        "--summary-out",
        type=Path,
        help="Write one phylogenetic signal summary ledger as TSV or CSV.",
    )
    comparative_signal.add_argument(
        "--permutations-out",
        type=Path,
        help="Write one permutation ledger for the Blomberg-K test as TSV or CSV.",
    )
    comparative_signal.add_argument("--permutations", type=int, default=199)
    comparative_signal.add_argument("--seed", type=int, default=1)
    comparative_signal.add_argument(
        "--json", action="store_true", help="Emit the signal report as JSON."
    )
    _add_manifest_argument(comparative_signal)


def run_phylogenetic_signal_command(args: Any) -> int | None:
    if args.comparative_command != "signal":
        return None

    report = summarize_phylogenetic_signal(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        permutations=args.permutations,
        seed=args.seed,
    )
    if args.summary_out:
        write_phylogenetic_signal_summary_table(args.summary_out, report)
    if args.permutations_out:
        write_phylogenetic_signal_permutation_table(
            args.permutations_out,
            report,
        )
    outputs = _finalize_outputs(
        args,
        command="comparative",
        inputs=[args.tree, args.table],
    )
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
            metrics={
                "taxon_count": report.taxon_count,
                "blombergs_k": report.blombergs_k.k,
                "pagels_lambda": report.pagels_lambda.lambda_value,
                "signal_p_value": report.signal_test.p_value,
                "tree_is_ultrametric": report.input_audit.tree_is_ultrametric,
                "ultrametric_policy": report.input_audit.ultrametric_policy,
                "missing_value_policy": report.input_audit.missing_value_policy,
                "pruned_missing_value_taxon_count": len(
                    report.input_audit.pruned_missing_value_taxa
                ),
                "signal_seed": report.signal_test.seed,
                "signal_null_k_minimum": report.signal_test.null_distribution_minimum,
                "signal_null_k_mean": report.signal_test.null_distribution_mean,
                "signal_null_k_maximum": report.signal_test.null_distribution_maximum,
                "lambda_log_likelihood": report.pagels_lambda.log_likelihood,
                "lambda_likelihood_ratio_statistic": (
                    report.lambda_likelihood_ratio_statistic
                ),
                "lambda_likelihood_ratio_p_value": (
                    report.lambda_likelihood_ratio_p_value
                ),
                "lambda_optimizer_name": (
                    report.pagels_lambda.optimizer_diagnostics.optimizer_name
                ),
                "lambda_optimizer_function_evaluation_count": (
                    report.pagels_lambda.optimizer_diagnostics.function_evaluation_count
                ),
                "lambda_optimizer_hit_lower_boundary": (
                    report.pagels_lambda.optimizer_diagnostics.hit_lower_boundary
                ),
                "lambda_optimizer_hit_upper_boundary": (
                    report.pagels_lambda.optimizer_diagnostics.hit_upper_boundary
                ),
                "permutation_row_count": len(report.signal_test.permutation_rows),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
