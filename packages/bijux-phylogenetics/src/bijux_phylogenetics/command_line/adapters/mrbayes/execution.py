from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    prepare_mrbayes_analysis,
    run_mrbayes_posterior_inference,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_external_adapter_execution_arguments,
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_mrbayes_execution_commands(adapter_subparsers: Any) -> None:
    adapter_mrbayes_prepare = adapter_subparsers.add_parser(
        "mrbayes-prepare",
        help="Prepare a MrBayes NEXUS analysis from an aligned FASTA file.",
    )
    adapter_mrbayes_prepare.add_argument("input_path", type=Path)
    adapter_mrbayes_prepare.add_argument("--out", required=True, type=Path)
    adapter_mrbayes_prepare.add_argument("--partitions", type=Path)
    adapter_mrbayes_prepare.add_argument("--model", default="gtr")
    adapter_mrbayes_prepare.add_argument("--rates", default="gamma")
    adapter_mrbayes_prepare.add_argument("--ngen", type=int, default=10000)
    adapter_mrbayes_prepare.add_argument("--nchains", type=int, default=4)
    adapter_mrbayes_prepare.add_argument("--samplefreq", type=int, default=100)
    adapter_mrbayes_prepare.add_argument("--printfreq", type=int, default=100)
    adapter_mrbayes_prepare.add_argument("--burnin-fraction", type=float, default=0.25)
    adapter_mrbayes_prepare.add_argument(
        "--json", action="store_true", help="Emit the preparation report as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_prepare)

    adapter_mrbayes_run = adapter_subparsers.add_parser(
        "mrbayes-run",
        help="Run a prepared MrBayes posterior inference workflow.",
    )
    adapter_mrbayes_run.add_argument("input_path", type=Path)
    adapter_mrbayes_run.add_argument("--executable", type=str)
    adapter_mrbayes_run.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_mrbayes_run)
    _add_manifest_argument(adapter_mrbayes_run)


def run_mrbayes_execution_command(args: Any) -> int | None:
    if args.adapter_command == "mrbayes-prepare":
        report = prepare_mrbayes_analysis(
            args.input_path,
            args.out,
            partition_path=args.partitions,
            model=args.model,
            rates=args.rates,
            ngen=args.ngen,
            nchains=args.nchains,
            samplefreq=args.samplefreq,
            printfreq=args.printfreq,
            burnin_fraction=args.burnin_fraction,
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
                    "partitioned": report.partition_path is not None,
                    "partition_count": report.partition_count,
                    "partition_warning_count": len(report.partition_warnings),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-run":
        report = run_mrbayes_posterior_inference(
            args.input_path,
            executable=args.executable or "mb",
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
                    "resumed": report.resumed,
                    "timeout_seconds": report.run.timeout_seconds,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
