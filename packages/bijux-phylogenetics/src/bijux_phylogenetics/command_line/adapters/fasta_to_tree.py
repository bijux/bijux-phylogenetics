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
    list_mafft_alignment_modes,
    list_trimal_trimming_modes,
    run_fasta_to_tree_workflow,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import method_tier_metrics
from bijux_phylogenetics.runtime.results import build_command_result


def add_adapter_fasta_to_tree_commands(adapter_subparsers: Any) -> None:
    adapter_fasta_to_tree = adapter_subparsers.add_parser(
        "fasta-to-tree", help="Run alignment-to-tree inference from raw FASTA."
    )
    adapter_fasta_to_tree.add_argument("input_path", type=Path)
    adapter_fasta_to_tree.add_argument("--out-dir", required=True, type=Path)
    adapter_fasta_to_tree.add_argument("--prefix")
    adapter_fasta_to_tree.add_argument(
        "--sequence-type", choices=("dna", "rna", "protein", "unknown")
    )
    adapter_fasta_to_tree.add_argument("--mafft-executable", type=str)
    adapter_fasta_to_tree.add_argument(
        "--alignment-mode",
        choices=list_mafft_alignment_modes(),
        default="auto",
        help="Select the named MAFFT alignment strategy for the raw-input workflow.",
    )
    adapter_fasta_to_tree.add_argument("--trimal-executable", type=str)
    adapter_fasta_to_tree.add_argument(
        "--trimming-mode",
        choices=list_trimal_trimming_modes(),
        default="gap-threshold",
        help="Select the named trimAl trimming strategy for the aligned workflow.",
    )
    adapter_fasta_to_tree.add_argument("--iqtree-executable", type=str)
    adapter_fasta_to_tree.add_argument(
        "--iqtree-seed",
        type=int,
        default=1,
        help="Set the IQ-TREE random seed for deterministic model selection and support estimation.",
    )
    adapter_fasta_to_tree.add_argument(
        "--iqtree-threads",
        type=int,
        default=1,
        help="Set the IQ-TREE thread count used across model selection and inference.",
    )
    adapter_fasta_to_tree.add_argument("--trim-gap-threshold", type=float, default=0.1)
    adapter_fasta_to_tree.add_argument("--bootstrap-replicates", type=int, default=1000)
    adapter_fasta_to_tree.add_argument(
        "--normalize-identifiers",
        action="store_true",
        help="Normalize FASTA identifiers before alignment and resolve any collisions.",
    )
    adapter_fasta_to_tree.add_argument(
        "--remove-invalid-records",
        action="store_true",
        help="Remove empty or illegal sequence records before alignment.",
    )
    adapter_fasta_to_tree.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_fasta_to_tree)
    _add_manifest_argument(adapter_fasta_to_tree)


def run_adapter_fasta_to_tree_command(args: Any) -> int | None:
    if args.adapter_command != "fasta-to-tree":
        return None

    report = run_fasta_to_tree_workflow(
        args.input_path,
        out_dir=args.out_dir,
        prefix=args.prefix,
        sequence_type=args.sequence_type,
        mafft_executable=args.mafft_executable or "mafft",
        alignment_mode=args.alignment_mode,
        trimal_executable=args.trimal_executable or "trimal",
        trimming_mode=args.trimming_mode,
        iqtree_executable=args.iqtree_executable or "iqtree2",
        iqtree_seed=args.iqtree_seed,
        iqtree_threads=args.iqtree_threads,
        trim_gap_threshold=args.trim_gap_threshold,
        bootstrap_replicates=args.bootstrap_replicates,
        normalize_identifiers=args.normalize_identifiers,
        remove_invalid_records=args.remove_invalid_records,
        resume=args.resume,
        timeout_seconds=args.timeout_seconds,
        incomplete_run_policy=args.incomplete_run_policy,
    )
    outputs = _finalize_outputs(
        args,
        command="adapter",
        inputs=[args.input_path],
        outputs=[report.engine_artifact_dir, *report.output_paths.values()],
    )
    _print_result(
        build_command_result(
            command="adapter",
            inputs=[args.input_path],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "alignment_mode": args.alignment_mode,
                "trimming_mode": args.trimming_mode,
                "iqtree_seed": args.iqtree_seed,
                "iqtree_threads": args.iqtree_threads,
                "bootstrap_replicates": args.bootstrap_replicates,
                "retained_site_count": (
                    None
                    if report.trimming_workflow.trimming_summary is None
                    else report.trimming_workflow.trimming_summary.retained_site_count
                ),
                "removed_site_count": (
                    None
                    if report.trimming_workflow.trimming_summary is None
                    else report.trimming_workflow.trimming_summary.removed_site_count
                ),
                "selected_model": report.selected_model,
                "sequence_type": report.sequence_type,
                "sequence_type_confidence": (
                    report.input_validation.sequence_type_report.confidence
                    if report.repaired_input_validation is None
                    else report.repaired_input_validation.sequence_type_report.confidence
                ),
                "normalized_identifier_count": (
                    0
                    if report.input_repair is None
                    else len(report.input_repair.normalized_identifiers)
                ),
                "removed_record_count": (
                    0
                    if report.input_repair is None
                    else len(report.input_repair.removed_records)
                ),
                "resumed": any(
                    workflow.resumed
                    for workflow in (
                        report.alignment_workflow,
                        report.trimming_workflow,
                        report.model_selection_workflow,
                        report.maximum_likelihood_workflow,
                        report.bootstrap_workflow,
                    )
                ),
                "timeout_seconds": args.timeout_seconds,
                **method_tier_metrics(report.method_tier),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
