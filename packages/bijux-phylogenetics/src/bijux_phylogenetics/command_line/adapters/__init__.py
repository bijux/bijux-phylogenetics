from __future__ import annotations

from pathlib import Path
from typing import Any

from .bayesian import (
    add_bayesian_adapter_commands,
    run_bayesian_adapter_command,
)
from .beast import (
    add_beast_adapter_commands,
    run_beast_adapter_command,
)
from .inference import (
    add_inference_adapter_commands,
    run_inference_adapter_command,
)
from .inspection import (
    add_adapter_inspection_commands,
    run_adapter_inspection_command,
)
from .mrbayes import (
    add_mrbayes_adapter_commands,
    run_mrbayes_adapter_command,
)
from .reporting import (
    add_adapter_reporting_commands,
    run_adapter_reporting_command,
)
from bijux_phylogenetics.command_line.arguments import (
    _add_external_adapter_execution_arguments,
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.engines import (
    list_mafft_alignment_modes,
    list_trimal_trimming_modes,
    run_alignment_trimming,
    run_bootstrap_consensus_tree,
    run_bootstrap_support_estimation,
    run_codon_aware_multiple_sequence_alignment,
    run_fasta_to_tree_workflow,
    run_maximum_likelihood_tree_inference,
    run_model_selection,
    run_multiple_sequence_alignment,
    run_sh_alrt_support_estimation,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import method_tier_metrics
from bijux_phylogenetics.runtime.errors import EngineUnavailableError
from bijux_phylogenetics.runtime.results import build_command_result


def add_adapter_commands(subparsers: Any) -> None:
    adapter = subparsers.add_parser(
        get_command_spec("adapter").name, help=get_command_spec("adapter").summary
    )
    adapter_subparsers = adapter.add_subparsers(dest="adapter_command", required=True)
    add_adapter_inspection_commands(adapter_subparsers)
    add_adapter_reporting_commands(adapter_subparsers)

    adapter_align = adapter_subparsers.add_parser(
        "align", help="Run multiple-sequence alignment on unaligned FASTA."
    )
    adapter_align.add_argument("input_path", type=Path)
    adapter_align.add_argument("--out", required=True, type=Path)
    adapter_align.add_argument("--executable", type=str)
    adapter_align.add_argument(
        "--mode",
        choices=list_mafft_alignment_modes(),
        default="auto",
        help="Select the named MAFFT alignment strategy.",
    )
    adapter_align.add_argument(
        "--codon-aware",
        action="store_true",
        help="Translate accepted coding nucleotide sequences to an amino-acid guide, align that guide, then back-translate codon triplets.",
    )
    adapter_align.add_argument(
        "--sequence-type",
        choices=("dna", "rna"),
        help="Declare the coding nucleotide type for codon-aware alignment when explicit forcing is needed.",
    )
    adapter_align.add_argument(
        "--genetic-code",
        default="1",
        help="Use an NCBI genetic code id or codon-table name for codon-aware validation and translation.",
    )
    adapter_align.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_align)
    _add_manifest_argument(adapter_align)

    adapter_trim = adapter_subparsers.add_parser(
        "trim", help="Run external alignment trimming."
    )
    adapter_trim.add_argument("input_path", type=Path)
    adapter_trim.add_argument("--out", required=True, type=Path)
    adapter_trim.add_argument(
        "--mode",
        choices=list_trimal_trimming_modes(),
        default="gap-threshold",
        help="Select the named trimAl trimming strategy.",
    )
    adapter_trim.add_argument("--gap-threshold", type=float, default=0.1)
    adapter_trim.add_argument("--executable", type=str)
    adapter_trim.add_argument(
        "--json", action="store_true", help="Emit the workflow report as JSON."
    )
    _add_external_adapter_execution_arguments(adapter_trim)
    _add_manifest_argument(adapter_trim)

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

    add_inference_adapter_commands(adapter_subparsers)
    add_mrbayes_adapter_commands(adapter_subparsers)
    add_beast_adapter_commands(adapter_subparsers)
    add_bayesian_adapter_commands(adapter_subparsers)


def run_adapter_command(args: Any) -> int | None:
    if args.command != "adapter":
        return None
    inspection_exit_code = run_adapter_inspection_command(args)
    if inspection_exit_code is not None:
        return inspection_exit_code
    reporting_exit_code = run_adapter_reporting_command(args)
    if reporting_exit_code is not None:
        return reporting_exit_code
    if args.adapter_command == "align":
        if args.codon_aware:
            report = run_codon_aware_multiple_sequence_alignment(
                args.input_path,
                args.out,
                executable=args.executable or "mafft",
                mode=args.mode,
                sequence_type=args.sequence_type,
                genetic_code=args.genetic_code,
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
                    warnings=report.warnings,
                    metrics={
                        "mode": args.mode,
                        "codon_aware": True,
                        "sequence_type": report.sequence_type,
                        "genetic_code_id": report.genetic_code_id,
                        "accepted_sequence_count": report.accepted_sequence_count,
                        "excluded_sequence_count": len(report.excluded_sequences),
                        "invalid_codon_sequence_count": report.invalid_codon_sequence_count,
                        "terminal_stop_sequence_count": report.terminal_stop_sequence_count,
                        "resumed": report.resumed,
                        "timeout_seconds": report.run.timeout_seconds,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        report = run_multiple_sequence_alignment(
            args.input_path,
            args.out,
            executable=args.executable or "mafft",
            mode=args.mode,
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
                    "mode": args.mode,
                    "codon_aware": False,
                    "warning_count": len(report.run.warning_lines),
                    "resumed": report.resumed,
                    "timeout_seconds": report.run.timeout_seconds,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.adapter_command == "trim":
        report = run_alignment_trimming(
            args.input_path,
            args.out,
            executable=args.executable or "trimal",
            mode=args.mode,
            gap_threshold=args.gap_threshold,
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
                    "mode": args.mode,
                    "warning_count": len(report.run.warning_lines),
                    "retained_site_count": (
                        None
                        if report.trimming_summary is None
                        else report.trimming_summary.retained_site_count
                    ),
                    "removed_site_count": (
                        None
                        if report.trimming_summary is None
                        else report.trimming_summary.removed_site_count
                    ),
                    "input_gap_percentage": (
                        None
                        if report.trimming_summary is None
                        else report.trimming_summary.input_gap_percentage
                    ),
                    "trimmed_gap_percentage": (
                        None
                        if report.trimming_summary is None
                        else report.trimming_summary.trimmed_gap_percentage
                    ),
                    "resumed": report.resumed,
                    "timeout_seconds": report.run.timeout_seconds,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
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
        adapter_inputs = (
            [args.input_path] if args.partitions is None else [args.input_path, args.partitions]
        )
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
    if args.adapter_command == "infer-ml":
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
        adapter_inputs = (
            [args.input_path] if args.partitions is None else [args.input_path, args.partitions]
        )
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
        adapter_inputs = (
            [args.input_path] if args.partitions is None else [args.input_path, args.partitions]
        )
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
        adapter_inputs = (
            [args.input_path] if args.partitions is None else [args.input_path, args.partitions]
        )
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
    if args.adapter_command == "fasta-to-tree":
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
    if args.adapter_command == "consensus":
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
    mrbayes_exit_code = run_mrbayes_adapter_command(args)
    if mrbayes_exit_code is not None:
        return mrbayes_exit_code
    beast_exit_code = run_beast_adapter_command(args)
    if beast_exit_code is not None:
        return beast_exit_code
    bayesian_exit_code = run_bayesian_adapter_command(args)
    if bayesian_exit_code is not None:
        return bayesian_exit_code
    inference_exit_code = run_inference_adapter_command(args)
    if inference_exit_code is not None:
        return inference_exit_code
    raise EngineUnavailableError(f"unsupported adapter command: {args.adapter_command}")
