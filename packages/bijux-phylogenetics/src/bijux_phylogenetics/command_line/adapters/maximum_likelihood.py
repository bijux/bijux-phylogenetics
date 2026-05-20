from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_external_adapter_execution_arguments,
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.engines import (
    run_maximum_likelihood_tree_inference,
    run_model_selection,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_adapter_maximum_likelihood_commands(adapter_subparsers: Any) -> None:
    adapter_model = adapter_subparsers.add_parser(
        "model-select", help="Run external sequence-model selection."
    )
    adapter_model.add_argument("input_path", type=Path)
    adapter_model.add_argument("--out-dir", required=True, type=Path)
    adapter_model.add_argument("--prefix", default="model-selection")
    adapter_model.add_argument(
        "--partitions",
        type=Path,
        help="Validate and apply a partition scheme for partitioned model selection.",
    )
    adapter_model.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_model.add_argument("--executable", type=str)
    adapter_model.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_model)
    _add_manifest_argument(adapter_model)

    adapter_ml = adapter_subparsers.add_parser(
        "infer-ml", help="Run maximum-likelihood tree inference."
    )
    adapter_ml.add_argument("input_path", type=Path)
    adapter_ml.add_argument("--out-dir", required=True, type=Path)
    adapter_ml.add_argument("--model", required=True)
    adapter_ml.add_argument("--prefix", default="maximum-likelihood")
    adapter_ml.add_argument(
        "--partitions",
        type=Path,
        help="Validate and apply a partition scheme for partitioned maximum-likelihood inference.",
    )
    adapter_ml.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_ml.add_argument("--executable", type=str)
    adapter_ml.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_ml)
    _add_manifest_argument(adapter_ml)


def run_adapter_maximum_likelihood_command(args: Any) -> int | None:
    if args.adapter_command == "model-select":
        report = run_model_selection(
            args.input_path,
            out_dir=args.out_dir,
            prefix=args.prefix,
            executable=args.executable or "iqtree2",
            sequence_type=args.sequence_type,
            partition_path=args.partitions,
            resume=args.resume,
            timeout_seconds=args.timeout_seconds,
            incomplete_run_policy=args.incomplete_run_policy,
        )
        adapter_inputs = _adapter_inputs(args)
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=adapter_inputs,
            outputs=[*report.output_paths.values(), report.manifest_path],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=adapter_inputs,
                outputs=outputs,
                warnings=report.run.warning_lines,
                metrics={
                    "selected_model": report.selected_model,
                    "selected_criterion": (
                        None
                        if report.model_selection_summary is None
                        else report.model_selection_summary.selected_criterion
                    ),
                    "candidate_model_count": (
                        0
                        if report.model_selection_summary is None
                        else report.model_selection_summary.candidate_count
                    ),
                    "best_model_aic": (
                        None
                        if report.model_selection_summary is None
                        else report.model_selection_summary.best_model_aic
                    ),
                    "best_model_aicc": (
                        None
                        if report.model_selection_summary is None
                        else report.model_selection_summary.best_model_aicc
                    ),
                    "best_model_bic": (
                        None
                        if report.model_selection_summary is None
                        else report.model_selection_summary.best_model_bic
                    ),
                    "log_likelihood": report.log_likelihood,
                    "partitioned": args.partitions is not None,
                    "resumed": report.resumed,
                    "timeout_seconds": report.run.timeout_seconds,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.adapter_command != "infer-ml":
        return None

    report = run_maximum_likelihood_tree_inference(
        args.input_path,
        out_dir=args.out_dir,
        model=args.model,
        prefix=args.prefix,
        executable=args.executable or "iqtree2",
        sequence_type=args.sequence_type,
        partition_path=args.partitions,
        resume=args.resume,
        timeout_seconds=args.timeout_seconds,
        incomplete_run_policy=args.incomplete_run_policy,
    )
    adapter_inputs = _adapter_inputs(args)
    outputs = _finalize_outputs(
        args,
        command="adapter",
        inputs=adapter_inputs,
        outputs=[*report.output_paths.values(), report.manifest_path],
    )
    _print_result(
        build_command_result(
            command="adapter",
            inputs=adapter_inputs,
            outputs=outputs,
            warnings=report.run.warning_lines,
            metrics={
                "selected_model": report.selected_model,
                "selected_criterion": (
                    None
                    if report.model_selection_summary is None
                    else report.model_selection_summary.selected_criterion
                ),
                "candidate_model_count": (
                    0
                    if report.model_selection_summary is None
                    else report.model_selection_summary.candidate_count
                ),
                "log_likelihood": report.log_likelihood,
                "support_value_count": (
                    0
                    if report.iqtree_summary is None
                    else report.iqtree_summary.support_value_count
                ),
                "partitioned": args.partitions is not None,
                "resumed": report.resumed,
                "timeout_seconds": report.run.timeout_seconds,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _adapter_inputs(args: Any) -> list[Path]:
    if args.partitions is None:
        return [args.input_path]
    return [args.input_path, args.partitions]
