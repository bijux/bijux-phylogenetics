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
    compare_fast_and_ml_trees,
    run_fast_tree_inference,
    run_tree_inference_comparison,
)
from bijux_phylogenetics.engines.inference import (
    run_inference_reproducibility_check,
    run_large_alignment_inference,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_inference_adapter_commands(adapter_subparsers: Any) -> None:
    adapter_fast = adapter_subparsers.add_parser(
        "infer-fast", help="Run fast approximate tree inference."
    )
    adapter_fast.add_argument("input_path", type=Path)
    adapter_fast.add_argument("--out", required=True, type=Path)
    adapter_fast.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_fast.add_argument("--executable", type=str)
    adapter_fast.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_fast)
    _add_manifest_argument(adapter_fast)

    adapter_large = adapter_subparsers.add_parser(
        "infer-large",
        help="Run large-alignment FastTree inference with streamed preflight and resource reporting.",
    )
    adapter_large.add_argument("input_path", type=Path)
    adapter_large.add_argument("--out-dir", required=True, type=Path)
    adapter_large.add_argument("--prefix", default="large-alignment-inference")
    adapter_large.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_large.add_argument("--executable", type=str)
    adapter_large.add_argument("--resume", action="store_true")
    adapter_large.add_argument(
        "--timeout-seconds",
        type=float,
        help="Stop the FastTree inference step if it exceeds this wall-clock budget.",
    )
    adapter_large.add_argument(
        "--incomplete-run-policy",
        choices=("reject", "clean"),
        default="reject",
        help="Reject or clean incomplete FastTree outputs before starting a new large-alignment run.",
    )
    adapter_large.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_manifest_argument(adapter_large)

    adapter_compare = adapter_subparsers.add_parser(
        "compare", help="Compare fast approximate and ML trees."
    )
    adapter_compare.add_argument("--fast-tree", required=True, type=Path)
    adapter_compare.add_argument("--ml-tree", required=True, type=Path)
    adapter_compare.add_argument("--out", required=True, type=Path)
    adapter_compare.add_argument(
        "--json", action="store_true", help="Emit the comparison report as JSON."
    )
    _add_manifest_argument(adapter_compare)

    adapter_compare_engines = adapter_subparsers.add_parser(
        "compare-engines",
        help="Run IQ-TREE and FastTree on one alignment and compare the inferred trees.",
    )
    adapter_compare_engines.add_argument("input_path", type=Path)
    adapter_compare_engines.add_argument("--out-dir", required=True, type=Path)
    adapter_compare_engines.add_argument("--prefix", default="engine-comparison")
    adapter_compare_engines.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_compare_engines.add_argument("--iqtree-executable", type=str)
    adapter_compare_engines.add_argument("--fasttree-executable", type=str)
    adapter_compare_engines.add_argument(
        "--iqtree-seed",
        type=int,
        default=1,
        help="Set the IQ-TREE random seed for deterministic comparison runs.",
    )
    adapter_compare_engines.add_argument(
        "--iqtree-threads",
        type=int,
        default=1,
        help="Set the IQ-TREE thread count used during the comparison run.",
    )
    adapter_compare_engines.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=1000,
        help="Set the ultrafast bootstrap replicate count used for the IQ-TREE support workflow.",
    )
    adapter_compare_engines.add_argument(
        "--json",
        action="store_true",
        help="Emit the comparison workflow report as JSON.",
    )
    _add_external_adapter_execution_arguments(adapter_compare_engines)
    _add_manifest_argument(adapter_compare_engines)

    adapter_reproducibility = adapter_subparsers.add_parser(
        "reproducibility",
        help="Rerun supported IQ-TREE inference and classify deterministic versus unstable outputs.",
    )
    adapter_reproducibility.add_argument("input_path", type=Path)
    adapter_reproducibility.add_argument("--out-dir", required=True, type=Path)
    adapter_reproducibility.add_argument(
        "--prefix", default="inference-reproducibility"
    )
    adapter_reproducibility.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_reproducibility.add_argument("--iqtree-executable", type=str)
    adapter_reproducibility.add_argument(
        "--iqtree-seed",
        type=int,
        default=1,
        help="Set the IQ-TREE random seed used for every rerun.",
    )
    adapter_reproducibility.add_argument(
        "--iqtree-threads",
        type=int,
        default=1,
        help="Set the IQ-TREE thread count used for every rerun.",
    )
    adapter_reproducibility.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=1000,
        help="Set the ultrafast bootstrap replicate count used for every rerun.",
    )
    adapter_reproducibility.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Set how many repeated supported-inference runs to compare.",
    )
    adapter_reproducibility.add_argument(
        "--json",
        action="store_true",
        help="Emit the reproducibility workflow report as JSON.",
    )
    _add_manifest_argument(adapter_reproducibility)


def run_inference_adapter_command(args: Any) -> int | None:
    if args.adapter_command == "infer-fast":
        report = run_fast_tree_inference(
            args.input_path,
            args.out,
            executable=args.executable or "FastTree",
            sequence_type=args.sequence_type,
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
                    "approximate_method": (
                        None
                        if report.fasttree_support_summary is None
                        else report.fasttree_support_summary.approximate_method
                    ),
                    "support_label_kind": (
                        None
                        if report.fasttree_support_summary is None
                        else report.fasttree_support_summary.support_label_kind
                    ),
                    "support_scale": (
                        None
                        if report.fasttree_support_summary is None
                        else report.fasttree_support_summary.support_scale
                    ),
                    "annotated_node_count": (
                        0
                        if report.fasttree_support_summary is None
                        else report.fasttree_support_summary.annotated_node_count
                    ),
                    "minimum_local_support": (
                        None
                        if report.fasttree_support_summary is None
                        else report.fasttree_support_summary.minimum_local_support
                    ),
                    "maximum_local_support": (
                        None
                        if report.fasttree_support_summary is None
                        else report.fasttree_support_summary.maximum_local_support
                    ),
                    "weakly_supported_clade_count": (
                        0
                        if report.fasttree_support_summary is None
                        else report.fasttree_support_summary.weakly_supported_clade_count
                    ),
                    "resumed": report.resumed,
                    "timeout_seconds": report.run.timeout_seconds,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "infer-large":
        report = run_large_alignment_inference(
            args.input_path,
            out_dir=args.out_dir,
            prefix=args.prefix,
            sequence_type=args.sequence_type,
            executable=args.executable or "FastTree",
            timeout_seconds=args.timeout_seconds,
            resume=args.resume,
            incomplete_run_policy=args.incomplete_run_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=[*report.output_paths.values()],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "sequence_count": report.input_summary.sequence_count,
                    "alignment_length": report.input_summary.alignment_length,
                    "total_site_cells": report.input_summary.total_site_cells,
                    "sequence_type": report.sequence_type,
                    "resumed": report.resumed,
                    "timeout_seconds": report.timeout_seconds,
                    "peak_memory_bytes": max(
                        (
                            row.peak_memory_bytes
                            for row in report.resource_rows
                            if row.peak_memory_bytes is not None
                        ),
                        default=None,
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "compare":
        report = compare_fast_and_ml_trees(
            args.fast_tree, args.ml_tree, out_path=args.out
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.fast_tree, args.ml_tree],
            outputs=[report.comparison_report.output_path],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.fast_tree, args.ml_tree],
                outputs=outputs,
                metrics={
                    "shared_taxa": len(report.comparison_report.topology.shared_taxa),
                    "robinson_foulds_distance": (
                        report.comparison_report.topology.robinson_foulds_distance
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "compare-engines":
        report = run_tree_inference_comparison(
            args.input_path,
            out_dir=args.out_dir,
            prefix=args.prefix,
            sequence_type=args.sequence_type,
            iqtree_executable=args.iqtree_executable or "iqtree2",
            fasttree_executable=args.fasttree_executable or "FastTree",
            iqtree_seed=args.iqtree_seed,
            iqtree_threads=args.iqtree_threads,
            bootstrap_replicates=args.bootstrap_replicates,
            resume=args.resume,
            timeout_seconds=args.timeout_seconds,
            incomplete_run_policy=args.incomplete_run_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=[*report.output_paths.values()],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "selected_model": report.selected_model,
                    "shared_taxa": len(report.engine_comparison.topology.shared_taxa),
                    "robinson_foulds_distance": (
                        report.engine_comparison.topology.robinson_foulds_distance
                    ),
                    "shared_clade_count": len(report.shared_clade_rows),
                    "conflicting_clade_count": len(report.conflicting_clade_rows),
                    "stable_clade_count": report.conclusion_summary.stable_clade_count,
                    "unstable_clade_count": report.conclusion_summary.unstable_clade_count,
                    "engine_specific_clade_count": (
                        report.conclusion_summary.engine_specific_clade_count
                    ),
                    "support_disagreement_count": sum(
                        1
                        for row in report.conflicting_clade_rows
                        if row.conflict_kind == "support_disagreement"
                    ),
                    "high_support_conflict_count": (
                        report.conclusion_summary.high_support_conflict_count
                    ),
                    "low_support_disagreement_count": (
                        report.conclusion_summary.low_support_disagreement_count
                    ),
                    "serious_conflict_count": (
                        report.conclusion_summary.serious_conflict_count
                    ),
                    "resumed": any(
                        workflow.resumed
                        for workflow in (
                            report.model_selection_workflow,
                            report.iqtree_support_workflow,
                            report.fasttree_workflow,
                        )
                    ),
                    "timeout_seconds": args.timeout_seconds,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "reproducibility":
        report = run_inference_reproducibility_check(
            args.input_path,
            out_dir=args.out_dir,
            prefix=args.prefix,
            sequence_type=args.sequence_type,
            executable=args.iqtree_executable or "iqtree2",
            repeats=args.repeats,
            bootstrap_replicates=args.bootstrap_replicates,
            seed=args.iqtree_seed,
            threads=args.iqtree_threads,
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=[*report.output_paths.values()],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "selected_model": report.selected_model,
                    "overall_status": report.overall_status,
                    "repeat_count": report.repeat_count,
                    "unstable_comparison_count": sum(
                        1
                        for row in report.comparison_rows
                        if row.classification == "unstable"
                    ),
                    "equivalent_comparison_count": sum(
                        1
                        for row in report.comparison_rows
                        if row.classification == "equivalent"
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
