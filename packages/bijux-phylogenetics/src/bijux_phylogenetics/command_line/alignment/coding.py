from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.fasta import write_fasta_alignment
from bijux_phylogenetics.io.fasta.coding import (
    inspect_coding_alignment,
    translate_coding_alignment,
    write_translation_codon_validation_table,
    write_translation_excluded_sequence_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_alignment_coding_commands(alignment_subparsers: Any) -> None:
    alignment_coding = alignment_subparsers.add_parser(
        "coding",
        help="Inspect a nucleotide coding alignment for frameshift-like lengths and stop codons.",
    )
    alignment_coding.add_argument("alignment", type=Path)
    alignment_coding.add_argument(
        "--genetic-code",
        default="1",
        help="Use an NCBI genetic code id or codon-table name for coding diagnostics.",
    )
    alignment_coding.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_coding)

    alignment_translate = alignment_subparsers.add_parser(
        "translate",
        help="Translate a nucleotide coding alignment to amino acids.",
    )
    alignment_translate.add_argument("alignment", type=Path)
    alignment_translate.add_argument("--out", required=True, type=Path)
    alignment_translate.add_argument(
        "--codon-validation-out",
        type=Path,
        help="Write a codon-level translation validation ledger as TSV.",
    )
    alignment_translate.add_argument(
        "--excluded-sequences-out",
        type=Path,
        help="Write any translation exclusions as TSV.",
    )
    alignment_translate.add_argument(
        "--genetic-code",
        default="1",
        help="Use an NCBI genetic code id or codon-table name for coding translation.",
    )
    alignment_translate.add_argument(
        "--json", action="store_true", help="Emit the translation report as JSON."
    )
    _add_manifest_argument(alignment_translate)


def run_alignment_coding_command(args: Any) -> int | None:
    if args.alignment_command == "coding":
        report = inspect_coding_alignment(
            args.alignment,
            genetic_code=args.genetic_code,
        )
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "genetic_code_id": report.genetic_code_id,
                    "frameshift_like_sequence_count": len(
                        report.frameshift_like_sequences
                    ),
                    "invalid_codon_count": len(report.invalid_codons),
                    "stop_codon_count": len(report.stop_codons),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.alignment_command == "translate":
        records, report = translate_coding_alignment(
            args.alignment,
            genetic_code=args.genetic_code,
        )
        output_path = write_fasta_alignment(args.out, records)
        outputs_to_finalize = [output_path]
        if args.codon_validation_out is not None:
            outputs_to_finalize.append(
                write_translation_codon_validation_table(
                    args.codon_validation_out,
                    report,
                )
            )
        if args.excluded_sequences_out is not None:
            outputs_to_finalize.append(
                write_translation_excluded_sequence_table(
                    args.excluded_sequences_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
            outputs=outputs_to_finalize,
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "genetic_code_id": report.genetic_code_id,
                    "translated_sequence_count": report.translated_sequence_count,
                    "translated_alignment_length": report.translated_alignment_length,
                    "dropped_trailing_nucleotide_count": report.dropped_trailing_nucleotide_count,
                    "invalid_codon_count": report.invalid_codon_count,
                    "stop_codon_count": report.stop_codon_count,
                    "internal_stop_sequence_count": report.internal_stop_sequence_count,
                    "terminal_stop_sequence_count": report.terminal_stop_sequence_count,
                    "trailing_partial_codon_sequence_count": report.trailing_partial_codon_sequence_count,
                    "excluded_sequence_count": len(report.excluded_sequences),
                },
                data=report,
                warnings=report.warnings,
            ),
            json_output=args.json,
        )
        return 0

    return None
