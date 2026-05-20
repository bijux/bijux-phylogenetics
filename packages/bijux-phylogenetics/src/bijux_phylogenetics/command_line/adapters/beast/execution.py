from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    prepare_beast_time_tree_analysis,
    run_beast_posterior_inference,
    summarize_beast_analysis_xml,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_external_adapter_execution_arguments,
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_beast_execution_commands(adapter_subparsers: Any) -> None:
    adapter_beast_prepare = adapter_subparsers.add_parser(
        "beast-prepare",
        help="Prepare a BEAST-style time-tree XML analysis from aligned sequences and dating inputs.",
    )
    adapter_beast_prepare.add_argument("input_path", type=Path)
    adapter_beast_prepare.add_argument("--out", required=True, type=Path)
    adapter_beast_prepare.add_argument("--tree", type=Path)
    adapter_beast_prepare.add_argument("--calibrations", type=Path)
    adapter_beast_prepare.add_argument("--tip-dates", type=Path)
    adapter_beast_prepare.add_argument("--clock-model", default="strict")
    adapter_beast_prepare.add_argument("--tree-prior", default="yule")
    adapter_beast_prepare.add_argument("--chain-length", type=int, default=1000000)
    adapter_beast_prepare.add_argument("--log-every", type=int, default=1000)
    adapter_beast_prepare.add_argument(
        "--json", action="store_true", help="Emit the preparation report as JSON."
    )
    _add_manifest_argument(adapter_beast_prepare)

    adapter_beast_xml = adapter_subparsers.add_parser(
        "beast-xml",
        help="Summarize and validate one prepared BEAST analysis XML.",
    )
    adapter_beast_xml.add_argument("input_path", type=Path)
    adapter_beast_xml.add_argument(
        "--json", action="store_true", help="Emit the XML summary report as JSON."
    )
    _add_manifest_argument(adapter_beast_xml)

    adapter_beast_run = adapter_subparsers.add_parser(
        "beast-run",
        help="Run a prepared BEAST posterior inference workflow.",
    )
    adapter_beast_run.add_argument("input_path", type=Path)
    adapter_beast_run.add_argument("--executable", type=str)
    adapter_beast_run.add_argument("--threads", type=int, default=1)
    adapter_beast_run.add_argument("--seed", type=int, default=1)
    adapter_beast_run.add_argument(
        "--no-overwrite",
        action="store_true",
        help="Keep any existing posterior outputs instead of passing the BEAST overwrite flag.",
    )
    adapter_beast_run.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_beast_run)
    _add_manifest_argument(adapter_beast_run)


def run_beast_execution_command(args: Any) -> int | None:
    if args.adapter_command == "beast-prepare":
        report = prepare_beast_time_tree_analysis(
            args.input_path,
            args.out,
            tree_path=args.tree,
            calibration_path=args.calibrations,
            tip_dates_path=args.tip_dates,
            clock_model=args.clock_model,
            tree_prior=args.tree_prior,
            chain_length=args.chain_length,
            log_every=args.log_every,
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=[args.out],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "taxon_count": report.taxon_count,
                    "character_count": report.character_count,
                    "calibration_count": report.calibration_count,
                    "tip_date_count": report.tip_date_count,
                    "warning_count": report.warning_count,
                    "starting_tree_source": report.starting_tree_source,
                    "beast_data_type": report.beast_data_type,
                    "substitution_model": report.substitution_model,
                    "clock_model": report.clock_model,
                    "tree_prior": report.tree_prior,
                    "chain_length": report.chain_length,
                    "log_every": report.log_every,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-xml":
        report = summarize_beast_analysis_xml(args.input_path)
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=[],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "valid": report.valid,
                    "issue_count": len(report.issues),
                    "taxon_count": report.taxon_count,
                    "character_count": report.character_count,
                    "calibration_count": report.calibration_count,
                    "tip_date_count": report.tip_date_count,
                    "chain_length": report.chain_length,
                    "logger_count": report.logger_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-run":
        report = run_beast_posterior_inference(
            args.input_path,
            executable=args.executable or "beast",
            overwrite=not args.no_overwrite,
            threads=args.threads,
            seed=args.seed,
            resume=args.resume,
            timeout_seconds=args.timeout_seconds,
            incomplete_run_policy=args.incomplete_run_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=[*report.output_paths.values(), report.manifest_path],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                warnings=report.run.warning_lines,
                metrics={
                    "warning_count": len(report.run.warning_lines),
                    "threads": args.threads,
                    "seed": args.seed,
                    "overwrite": not args.no_overwrite,
                    "resumed": report.resumed,
                    "timeout_seconds": report.run.timeout_seconds,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
