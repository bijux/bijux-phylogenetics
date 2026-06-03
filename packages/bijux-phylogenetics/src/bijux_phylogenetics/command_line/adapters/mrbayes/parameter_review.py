from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    assess_mrbayes_burnin_sensitivity,
    assess_mrbayes_convergence,
    summarize_mrbayes_parameter_diagnostics,
    summarize_mrbayes_posterior_decomposition,
    write_mrbayes_burnin_sensitivity_slice_table,
    write_mrbayes_parameter_summary_table,
    write_mrbayes_posterior_decomposition_table,
)
from bijux_phylogenetics.bayesian.posterior_sets.burnin import (
    DEFAULT_BURNIN_FRACTIONS,
    write_burnin_clade_shift_table,
    write_burnin_parameter_shift_table,
)
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_mrbayes_parameter_review_commands(adapter_subparsers: Any) -> None:
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
        "--decomposition-out",
        type=Path,
        help="Write a TSV posterior, likelihood, and prior decomposition table for the retained trace samples.",
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


def run_mrbayes_parameter_review_command(args: Any) -> int | None:
    if args.adapter_command == "mrbayes-parameters":
        report = summarize_mrbayes_parameter_diagnostics(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
        )
        decomposition = None
        if args.decomposition_out is not None:
            decomposition = summarize_mrbayes_posterior_decomposition(
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
        if args.decomposition_out is not None:
            outputs.append(
                write_mrbayes_posterior_decomposition_table(
                    args.decomposition_out,
                    decomposition,
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
                    **(
                        {}
                        if decomposition is None
                        else {
                            "decomposition_verified": decomposition.verified,
                            "decomposition_maximum_absolute_delta": (
                                decomposition.maximum_absolute_delta
                            ),
                        }
                    ),
                },
                data=(
                    report
                    if decomposition is None
                    else {"summary": report, "decomposition": decomposition}
                ),
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
        inputs = [
            args.posterior_trees,
            *([args.traces] if args.traces is not None else []),
        ]
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

    return None
