from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.fasta import detect_invalid_alignment_characters
from bijux_phylogenetics.io.fasta.matrix import (
    compute_alignment_base_frequency_report,
    compute_alignment_segregating_site_report,
    write_alignment_base_frequency_table,
    write_alignment_segregating_site_table,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.runtime.results import build_command_result


def add_alignment_composition_commands(alignment_subparsers: Any) -> None:
    alignment_composition = alignment_subparsers.add_parser(
        "composition",
        help="Inspect inferred alphabet, composition, and GC content.",
    )
    alignment_composition.add_argument("alignment", type=Path)
    alignment_composition.add_argument(
        "--base-frequency-out",
        type=Path,
        help="Write ape-style alignment and per-sequence nucleotide state frequencies as TSV.",
    )
    alignment_composition.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_composition)

    alignment_segregating_sites = alignment_subparsers.add_parser(
        "segregating-sites",
        help="Report ape-style segregating alignment columns and their site indices.",
    )
    alignment_segregating_sites.add_argument("alignment", type=Path)
    alignment_segregating_sites.add_argument(
        "--site-table-out",
        type=Path,
        help="Write one TSV ledger of segregating alignment sites.",
    )
    alignment_segregating_sites.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_segregating_sites)

    alignment_invalid = alignment_subparsers.add_parser(
        "invalid",
        help="List alignment characters invalid for a declared alphabet.",
    )
    alignment_invalid.add_argument("alignment", type=Path)
    alignment_invalid.add_argument(
        "--alphabet",
        choices=("dna", "rna", "protein"),
        required=True,
    )
    alignment_invalid.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_invalid)


def run_alignment_composition_command(args: Any) -> int | None:
    if args.alignment_command == "composition":
        report = summarise_fasta(args.alignment)
        base_frequency_report = None
        warnings: list[str] = []
        try:
            base_frequency_report = compute_alignment_base_frequency_report(
                args.alignment
            )
        except InvalidAlignmentError:
            base_frequency_report = None
        else:
            warnings.extend(base_frequency_report.warnings)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        if args.base_frequency_out is not None:
            if base_frequency_report is None:
                raise InvalidAlignmentError(
                    "ape-style nucleotide base frequencies require a dna or rna alignment"
                )
            outputs.append(
                write_alignment_base_frequency_table(
                    args.base_frequency_out,
                    base_frequency_report,
                )
            )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=warnings,
                metrics={
                    "alphabet": report.inferred_alphabet,
                    "sequence_count": report.sequence_count,
                    "gc_sequence_count": len(report.per_sequence_gc_content),
                    "composition_outlier_count": len(report.composition_outliers),
                    "base_frequency_state_count": (
                        0
                        if base_frequency_report is None
                        else len(base_frequency_report.alignment_rows)
                    ),
                },
                data={
                    "path": report.path,
                    "inferred_alphabet": report.inferred_alphabet,
                    "nucleotide_composition": report.nucleotide_composition,
                    "amino_acid_composition": report.amino_acid_composition,
                    "per_sequence_gc_content": report.per_sequence_gc_content,
                    "whole_alignment_gc_content": report.whole_alignment_gc_content,
                    "composition_outliers": report.composition_outliers,
                    "alignment_state_frequencies": (
                        []
                        if base_frequency_report is None
                        else base_frequency_report.alignment_rows
                    ),
                    "per_sequence_state_frequencies": (
                        []
                        if base_frequency_report is None
                        else base_frequency_report.per_sequence_rows
                    ),
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "segregating-sites":
        report = compute_alignment_segregating_site_report(args.alignment)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        if args.site_table_out is not None:
            outputs.append(
                write_alignment_segregating_site_table(
                    args.site_table_out,
                    report,
                )
            )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "alphabet": report.inferred_alphabet,
                    "sequence_count": report.sequence_count,
                    "alignment_length": report.alignment_length,
                    "segregating_site_count": len(report.segregating_site_positions),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "invalid":
        report = detect_invalid_alignment_characters(
            args.alignment,
            alphabet=args.alphabet,
        )
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "invalid_character_count": len(report),
                    "alphabet": args.alphabet,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    return None
