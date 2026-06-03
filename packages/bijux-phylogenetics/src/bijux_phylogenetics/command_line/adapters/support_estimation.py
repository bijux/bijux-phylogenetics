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
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_sh_alrt_support_estimation,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_adapter_support_estimation_commands(adapter_subparsers: Any) -> None:
    adapter_bootstrap = adapter_subparsers.add_parser(
        "bootstrap", help="Run bootstrap support estimation."
    )
    adapter_bootstrap.add_argument("input_path", type=Path)
    adapter_bootstrap.add_argument("--out-dir", required=True, type=Path)
    adapter_bootstrap.add_argument("--model", required=True)
    adapter_bootstrap.add_argument("--replicates", type=int, default=1000)
    adapter_bootstrap.add_argument("--prefix", default="bootstrap-support")
    adapter_bootstrap.add_argument(
        "--partitions",
        type=Path,
        help="Validate and apply a partition scheme for partitioned bootstrap support estimation.",
    )
    adapter_bootstrap.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_bootstrap.add_argument("--executable", type=str)
    adapter_bootstrap.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_bootstrap)
    _add_manifest_argument(adapter_bootstrap)

    adapter_sh_alrt = adapter_subparsers.add_parser(
        "sh-alrt",
        help="Run combined sh-alrt and ultrafast bootstrap support estimation.",
    )
    adapter_sh_alrt.add_argument("input_path", type=Path)
    adapter_sh_alrt.add_argument("--out-dir", required=True, type=Path)
    adapter_sh_alrt.add_argument("--model", required=True)
    adapter_sh_alrt.add_argument("--alrt-replicates", type=int, default=1000)
    adapter_sh_alrt.add_argument("--bootstrap-replicates", type=int, default=1000)
    adapter_sh_alrt.add_argument("--prefix", default="sh-alrt-support")
    adapter_sh_alrt.add_argument(
        "--partitions",
        type=Path,
        help="Validate and apply a partition scheme for combined sh-alrt and ultrafast bootstrap support estimation.",
    )
    adapter_sh_alrt.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_sh_alrt.add_argument("--executable", type=str)
    adapter_sh_alrt.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_sh_alrt)
    _add_manifest_argument(adapter_sh_alrt)

    adapter_consensus = adapter_subparsers.add_parser(
        "consensus", help="Build a consensus tree from bootstrap trees."
    )
    adapter_consensus.add_argument("input_path", type=Path)
    adapter_consensus.add_argument("--out-dir", required=True, type=Path)
    adapter_consensus.add_argument("--prefix", default="bootstrap-consensus")
    adapter_consensus.add_argument("--minimum-support", type=float, default=0.5)
    adapter_consensus.add_argument("--executable", type=str)
    adapter_consensus.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_consensus)
    _add_manifest_argument(adapter_consensus)


def run_adapter_support_estimation_command(args: Any) -> int | None:
    if args.adapter_command == "bootstrap":
        report = run_bootstrap_support_estimation(
            args.input_path,
            out_dir=args.out_dir,
            model=args.model,
            replicates=args.replicates,
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
                    "bootstrap_replicates": args.replicates,
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
                    "minimum_support": (
                        None
                        if report.bootstrap_support_summary is None
                        else report.bootstrap_support_summary.minimum_support
                    ),
                    "maximum_support": (
                        None
                        if report.bootstrap_support_summary is None
                        else report.bootstrap_support_summary.maximum_support
                    ),
                    "weakly_supported_clade_count": (
                        0
                        if report.bootstrap_support_summary is None
                        else report.bootstrap_support_summary.weakly_supported_clade_count
                    ),
                    "weak_backbone_node_count": (
                        0
                        if report.weak_backbone_report is None
                        else report.weak_backbone_report.weak_backbone_node_count
                    ),
                    "support_histogram": (
                        {}
                        if report.bootstrap_support_summary is None
                        else report.bootstrap_support_summary.support_histogram
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
    if args.adapter_command == "sh-alrt":
        report = run_sh_alrt_support_estimation(
            args.input_path,
            out_dir=args.out_dir,
            model=args.model,
            sh_alrt_replicates=args.alrt_replicates,
            bootstrap_replicates=args.bootstrap_replicates,
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
                    "sh_alrt_replicates": args.alrt_replicates,
                    "bootstrap_replicates": args.bootstrap_replicates,
                    "selected_model": report.selected_model,
                    "log_likelihood": report.log_likelihood,
                    "support_value_count": (
                        0
                        if report.iqtree_summary is None
                        else report.iqtree_summary.support_value_count
                    ),
                    "sh_alrt_supported_node_count": (
                        0
                        if report.sh_alrt_support_summary is None
                        else report.sh_alrt_support_summary.annotated_node_count
                    ),
                    "conflicting_support_signal_count": (
                        0
                        if report.sh_alrt_support_summary is None
                        else report.sh_alrt_support_summary.conflicting_support_signal_count
                    ),
                    "minimum_sh_alrt_support": (
                        None
                        if report.sh_alrt_support_summary is None
                        else report.sh_alrt_support_summary.minimum_sh_alrt_support
                    ),
                    "maximum_sh_alrt_support": (
                        None
                        if report.sh_alrt_support_summary is None
                        else report.sh_alrt_support_summary.maximum_sh_alrt_support
                    ),
                    "minimum_ufboot_support": (
                        None
                        if report.sh_alrt_support_summary is None
                        else report.sh_alrt_support_summary.minimum_ufboot_support
                    ),
                    "maximum_ufboot_support": (
                        None
                        if report.sh_alrt_support_summary is None
                        else report.sh_alrt_support_summary.maximum_ufboot_support
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
    if args.adapter_command != "consensus":
        return None

    report = run_bootstrap_consensus_tree(
        args.input_path,
        out_dir=args.out_dir,
        prefix=args.prefix,
        executable=args.executable or "iqtree2",
        minimum_support=args.minimum_support,
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
                "minimum_support": args.minimum_support,
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
