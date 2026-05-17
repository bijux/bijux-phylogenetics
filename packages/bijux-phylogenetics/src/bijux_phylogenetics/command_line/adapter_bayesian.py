from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    build_bayesian_evidence_package,
    write_bayesian_methods_summary_text,
    write_supplementary_bayesian_diagnostics_table,
)
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_bayesian_adapter_commands(adapter_subparsers: Any) -> None:
    adapter_bayesian_evidence = adapter_subparsers.add_parser(
        "bayesian-evidence",
        help="Bundle Bayesian configs, trees, logs, diagnostics, and reports into one evidence package.",
    )
    adapter_bayesian_evidence.add_argument("--out-dir", required=True, type=Path)
    adapter_bayesian_evidence.add_argument(
        "--inputs", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--configs", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--trees", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--logs", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--diagnostics", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--reports", nargs="+", required=True, type=Path
    )
    adapter_bayesian_evidence.add_argument(
        "--json", action="store_true", help="Emit the evidence-package report as JSON."
    )
    _add_manifest_argument(adapter_bayesian_evidence)

    adapter_bayesian_table = adapter_subparsers.add_parser(
        "bayesian-diagnostics-table",
        help="Write a supplementary Bayesian diagnostics table from posterior logs.",
    )
    adapter_bayesian_table.add_argument("posterior_trees", type=Path)
    adapter_bayesian_table.add_argument("--log", required=True, type=Path)
    adapter_bayesian_table.add_argument("--additional-logs", nargs="*", type=Path)
    adapter_bayesian_table.add_argument("--out", required=True, type=Path)
    adapter_bayesian_table.add_argument(
        "--burnin-fractions", nargs="+", type=float, default=[0.1, 0.25, 0.5]
    )
    adapter_bayesian_table.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_bayesian_table.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_bayesian_table.add_argument(
        "--cross-chain-mean-shift-threshold", type=float, default=0.75
    )
    adapter_bayesian_table.add_argument(
        "--json", action="store_true", help="Emit the diagnostics-table result as JSON."
    )
    _add_manifest_argument(adapter_bayesian_table)

    adapter_bayesian_methods = adapter_subparsers.add_parser(
        "bayesian-methods",
        help="Write reviewer-facing Bayesian methods summary text.",
    )
    adapter_bayesian_methods.add_argument("posterior_trees", type=Path)
    adapter_bayesian_methods.add_argument("--log", required=True, type=Path)
    adapter_bayesian_methods.add_argument("--additional-logs", nargs="*", type=Path)
    adapter_bayesian_methods.add_argument("--analysis-xml", type=Path)
    adapter_bayesian_methods.add_argument("--out", required=True, type=Path)
    adapter_bayesian_methods.add_argument("--tree-prior", default="unspecified")
    adapter_bayesian_methods.add_argument("--clock-model", default="unspecified")
    adapter_bayesian_methods.add_argument("--calibration-path", type=Path)
    adapter_bayesian_methods.add_argument("--tip-dates-path", type=Path)
    adapter_bayesian_methods.add_argument(
        "--burnin-fractions", nargs="+", type=float, default=[0.1, 0.25, 0.5]
    )
    adapter_bayesian_methods.add_argument("--ess-threshold", type=float, default=200.0)
    adapter_bayesian_methods.add_argument(
        "--mean-shift-threshold", type=float, default=0.5
    )
    adapter_bayesian_methods.add_argument(
        "--cross-chain-mean-shift-threshold", type=float, default=0.75
    )
    adapter_bayesian_methods.add_argument(
        "--json", action="store_true", help="Emit the methods-summary result as JSON."
    )
    _add_manifest_argument(adapter_bayesian_methods)


def run_bayesian_adapter_command(args: Any) -> int | None:
    if not str(args.adapter_command).startswith("bayesian-"):
        return None

    if args.adapter_command == "bayesian-evidence":
        report = build_bayesian_evidence_package(
            bundle_root=args.out_dir,
            input_paths=args.inputs,
            config_paths=args.configs,
            tree_paths=args.trees,
            log_paths=args.logs,
            diagnostic_paths=args.diagnostics,
            report_paths=args.reports,
        )
        inputs = [
            *args.inputs,
            *args.configs,
            *args.trees,
            *args.logs,
            *args.diagnostics,
            *args.reports,
        ]
        outputs = _finalize_outputs(
            args, command="adapter", inputs=inputs, outputs=[args.out_dir]
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=inputs,
                outputs=outputs,
                metrics={
                    "file_count": report.file_count,
                    "valid": report.valid,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "bayesian-diagnostics-table":
        report = write_supplementary_bayesian_diagnostics_table(
            args.out,
            posterior_tree_path=args.posterior_trees,
            primary_log_path=args.log,
            additional_log_paths=args.additional_logs,
            burnin_fractions=tuple(args.burnin_fractions),
            ess_threshold=args.ess_threshold,
            mean_shift_threshold=args.mean_shift_threshold,
            cross_chain_mean_shift_threshold=args.cross_chain_mean_shift_threshold,
        )
        inputs = [args.posterior_trees, args.log, *(args.additional_logs or [])]
        outputs = _finalize_outputs(
            args, command="adapter", inputs=inputs, outputs=[args.out]
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=inputs,
                outputs=outputs,
                metrics={
                    "row_count": report.row_count,
                    "warning_count": report.warning_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "bayesian-methods":
        report = write_bayesian_methods_summary_text(
            args.out,
            posterior_tree_path=args.posterior_trees,
            primary_log_path=args.log,
            additional_log_paths=args.additional_logs,
            analysis_xml_path=args.analysis_xml,
            tree_prior=args.tree_prior,
            clock_model=args.clock_model,
            calibration_path=args.calibration_path,
            tip_dates_path=args.tip_dates_path,
            burnin_fractions=tuple(args.burnin_fractions),
            ess_threshold=args.ess_threshold,
            mean_shift_threshold=args.mean_shift_threshold,
            cross_chain_mean_shift_threshold=args.cross_chain_mean_shift_threshold,
        )
        inputs = [
            args.posterior_trees,
            args.log,
            *(args.additional_logs or []),
            *([args.analysis_xml] if args.analysis_xml is not None else []),
        ]
        outputs = _finalize_outputs(
            args, command="adapter", inputs=inputs, outputs=[args.out]
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=inputs,
                outputs=outputs,
                metrics={"warning_count": report.warning_count},
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
