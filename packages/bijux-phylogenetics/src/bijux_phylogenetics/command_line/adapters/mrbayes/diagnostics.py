from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    compute_mrbayes_effective_sample_sizes,
    parse_mrbayes_consensus_tree,
    parse_mrbayes_mcmc_diagnostics,
    parse_mrbayes_parameter_traces,
)
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_mrbayes_diagnostic_commands(adapter_subparsers: Any) -> None:
    adapter_mrbayes_traces = adapter_subparsers.add_parser(
        "mrbayes-traces",
        help="Parse a MrBayes parameter trace table.",
    )
    adapter_mrbayes_traces.add_argument("input_path", type=Path)
    adapter_mrbayes_traces.add_argument(
        "--json", action="store_true", help="Emit the trace report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_traces)

    adapter_mrbayes_mcmc = adapter_subparsers.add_parser(
        "mrbayes-mcmc",
        help="Parse a MrBayes MCMC diagnostics table.",
    )
    adapter_mrbayes_mcmc.add_argument("input_path", type=Path)
    adapter_mrbayes_mcmc.add_argument(
        "--json", action="store_true", help="Emit the MCMC diagnostics report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_mcmc)

    adapter_mrbayes_consensus = adapter_subparsers.add_parser(
        "mrbayes-consensus",
        help="Parse a MrBayes consensus tree with posterior-probability annotations.",
    )
    adapter_mrbayes_consensus.add_argument("input_path", type=Path)
    adapter_mrbayes_consensus.add_argument(
        "--json", action="store_true", help="Emit the consensus tree report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_consensus)

    adapter_mrbayes_ess = adapter_subparsers.add_parser(
        "mrbayes-ess",
        help="Compute effective sample sizes from a MrBayes trace table.",
    )
    adapter_mrbayes_ess.add_argument("input_path", type=Path)
    adapter_mrbayes_ess.add_argument(
        "--json", action="store_true", help="Emit the ESS report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_ess)


def run_mrbayes_diagnostic_command(args: Any) -> int | None:
    if args.adapter_command == "mrbayes-traces":
        report = parse_mrbayes_parameter_traces(args.input_path)
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "row_count": report.row_count,
                    "column_count": len(report.columns),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-mcmc":
        report = parse_mrbayes_mcmc_diagnostics(args.input_path)
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "row_count": report.row_count,
                    "column_count": len(report.columns),
                    "comment_count": len(report.comment_lines),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-consensus":
        tree, report = parse_mrbayes_consensus_tree(args.input_path)
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "tip_count": tree.tip_count,
                    "annotated_node_count": report.annotated_node_count,
                    "maximum_posterior_probability": (
                        report.maximum_posterior_probability
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-ess":
        report = compute_mrbayes_effective_sample_sizes(args.input_path)
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={"parameter_count": len(report.effective_sample_sizes)},
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
