from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    assess_mrbayes_burnin_sensitivity,
    assess_mrbayes_convergence,
    compute_mrbayes_effective_sample_sizes,
    parse_mrbayes_consensus_tree,
    parse_mrbayes_mcmc_diagnostics,
    parse_mrbayes_parameter_traces,
    render_bayesian_posterior_report,
    summarize_mrbayes_parameter_diagnostics,
    write_mrbayes_burnin_sensitivity_slice_table,
    write_mrbayes_parameter_summary_table,
)
from bijux_phylogenetics.bayesian.posterior_sets.burnin import (
    DEFAULT_BURNIN_FRACTIONS,
    write_burnin_clade_shift_table,
    write_burnin_parameter_shift_table,
)
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    method_tier_metrics,
    method_tier_warnings,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .execution import (
    add_mrbayes_execution_commands,
    run_mrbayes_execution_command,
)
from .posterior_trees import (
    add_mrbayes_posterior_tree_commands,
    run_mrbayes_posterior_tree_command,
)


def add_mrbayes_adapter_commands(adapter_subparsers: Any) -> None:
    add_mrbayes_execution_commands(adapter_subparsers)
    add_mrbayes_posterior_tree_commands(adapter_subparsers)
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

    adapter_mrbayes_parameters = adapter_subparsers.add_parser(
        "mrbayes-parameters",
        help="Summarize burn-in-aware posterior parameter diagnostics from a MrBayes trace table.",
    )
    adapter_mrbayes_parameters.add_argument("input_path", type=Path)
    adapter_mrbayes_parameters.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early samples before reporting posterior summaries.",
    )
    adapter_mrbayes_parameters.add_argument(
        "--summary-out",
        type=Path,
        help="Write a TSV parameter-summary table for the retained trace samples.",
    )
    adapter_mrbayes_parameters.add_argument(
        "--json", action="store_true", help="Emit the parameter diagnostics as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_parameters)

    adapter_mrbayes_burnin = adapter_subparsers.add_parser(
        "mrbayes-burnin-sensitivity",
        help="Compare MrBayes posterior summaries across multiple burn-in fractions.",
    )
    adapter_mrbayes_burnin.add_argument("posterior_trees", type=Path)
    adapter_mrbayes_burnin.add_argument("--traces", type=Path)
    adapter_mrbayes_burnin.add_argument(
        "--burnin-fractions",
        nargs="+",
        type=float,
        default=list(DEFAULT_BURNIN_FRACTIONS),
    )
    adapter_mrbayes_burnin.add_argument("--slice-out", type=Path)
    adapter_mrbayes_burnin.add_argument("--parameter-out", type=Path)
    adapter_mrbayes_burnin.add_argument("--clade-out", type=Path)
    adapter_mrbayes_burnin.add_argument(
        "--json",
        action="store_true",
        help="Emit the burn-in sensitivity report as JSON.",
    )
    _add_manifest_argument(adapter_mrbayes_burnin)

    adapter_mrbayes_convergence = adapter_subparsers.add_parser(
        "mrbayes-convergence",
        help="Assess MrBayes trace convergence from ESS and trace drift.",
    )
    adapter_mrbayes_convergence.add_argument("input_path", type=Path)
    adapter_mrbayes_convergence.add_argument(
        "--ess-threshold", type=float, default=200.0
    )
    adapter_mrbayes_convergence.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_mrbayes_convergence.add_argument(
        "--json", action="store_true", help="Emit the convergence report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_convergence)

    adapter_mrbayes_report = adapter_subparsers.add_parser(
        "mrbayes-report",
        help="Render an HTML Bayesian posterior report from posterior trees and traces.",
    )
    adapter_mrbayes_report.add_argument("posterior_trees", type=Path)
    adapter_mrbayes_report.add_argument("--traces", required=True, type=Path)
    adapter_mrbayes_report.add_argument("--out", required=True, type=Path)
    adapter_mrbayes_report.add_argument("--burnin-fraction", type=float, default=0.25)
    adapter_mrbayes_report.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_mrbayes_report.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_mrbayes_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_report)


def run_mrbayes_adapter_command(args: Any) -> int | None:
    if not str(args.adapter_command).startswith("mrbayes-"):
        return None

    execution_result = run_mrbayes_execution_command(args)
    if execution_result is not None:
        return execution_result

    posterior_tree_result = run_mrbayes_posterior_tree_command(args)
    if posterior_tree_result is not None:
        return posterior_tree_result

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

    if args.adapter_command == "mrbayes-parameters":
        report = summarize_mrbayes_parameter_diagnostics(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
        )
        outputs: list[Path | str] = []
        if args.summary_out is not None:
            outputs.append(
                write_mrbayes_parameter_summary_table(
                    args.summary_out,
                    report,
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
                    "burnin_fraction": report.burnin_fraction,
                    "kept_row_count": report.kept_row_count,
                    "parameter_count": len(report.parameter_summaries),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-burnin-sensitivity":
        report = assess_mrbayes_burnin_sensitivity(
            args.posterior_trees,
            trace_path=args.traces,
            burnin_fractions=tuple(args.burnin_fractions),
        )
        inputs = [args.posterior_trees, *([args.traces] if args.traces is not None else [])]
        outputs: list[Path | str] = []
        if args.slice_out is not None:
            outputs.append(
                write_mrbayes_burnin_sensitivity_slice_table(
                    args.slice_out,
                    report,
                )
            )
        if args.parameter_out is not None:
            outputs.append(
                write_burnin_parameter_shift_table(
                    args.parameter_out,
                    report.parameter_shifts,
                )
            )
        if args.clade_out is not None:
            outputs.append(
                write_burnin_clade_shift_table(
                    args.clade_out,
                    report.clade_shifts,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=inputs,
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=inputs,
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "slice_count": len(report.slices),
                    "parameter_shift_count": len(report.parameter_shifts),
                    "unstable_parameter_count": report.unstable_parameter_count,
                    "clade_shift_count": len(report.clade_shifts),
                    "unstable_clade_count": report.unstable_clade_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-convergence":
        report = assess_mrbayes_convergence(
            args.input_path,
            ess_threshold=args.ess_threshold,
            mean_shift_threshold=args.mean_shift_threshold,
        )
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                warnings=[warning["message"] for warning in report.warnings],
                metrics={
                    "warning_count": len(report.warnings),
                    "converged": report.converged,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-report":
        report = render_bayesian_posterior_report(
            posterior_tree_path=args.posterior_trees,
            trace_path=args.traces,
            out_path=args.out,
            burnin_fraction=args.burnin_fraction,
            ess_threshold=args.ess_threshold,
            mean_shift_threshold=args.mean_shift_threshold,
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.posterior_trees, args.traces],
            outputs=[report.output_path],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.posterior_trees, args.traces],
                outputs=outputs,
                warnings=method_tier_warnings(report.method_tier),
                metrics={
                    "kept_tree_count": report.kept_tree_count,
                    "warning_count": report.warning_count
                    + len(method_tier_warnings(report.method_tier)),
                    **method_tier_metrics(report.method_tier),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
