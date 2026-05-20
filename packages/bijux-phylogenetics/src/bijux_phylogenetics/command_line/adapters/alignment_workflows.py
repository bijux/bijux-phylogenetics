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
    run_alignment_trimming,
    run_codon_aware_multiple_sequence_alignment,
    run_multiple_sequence_alignment,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_adapter_alignment_workflow_commands(adapter_subparsers: Any) -> None:
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


def run_adapter_alignment_workflow_command(args: Any) -> int | None:
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
    if args.adapter_command != "trim":
        return None

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
