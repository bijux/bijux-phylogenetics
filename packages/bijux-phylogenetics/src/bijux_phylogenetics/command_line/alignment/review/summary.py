from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.fasta import infer_alignment_alphabet, load_fasta_alignment
from bijux_phylogenetics.io.fasta.cleaning import list_alignment_filter_profiles
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.runtime.results import build_command_result


def add_alignment_summary_commands(alignment_subparsers: Any) -> None:
    alignment_alphabet = alignment_subparsers.add_parser(
        "alphabet",
        help="Infer the alignment sequence alphabet.",
    )
    alignment_alphabet.add_argument("alignment", type=Path)
    alignment_alphabet.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_alphabet)

    alignment_profiles = alignment_subparsers.add_parser(
        "profiles",
        help="List the supported named alignment-cleaning profiles.",
    )
    alignment_profiles.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_profiles)

    alignment_gc = alignment_subparsers.add_parser(
        "gc",
        help="Report per-sequence and whole-alignment GC content.",
    )
    alignment_gc.add_argument("alignment", type=Path)
    alignment_gc.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_gc)

    alignment_inspect = alignment_subparsers.add_parser(
        "inspect",
        help="Inspect an aligned FASTA file.",
    )
    alignment_inspect.add_argument("alignment", type=Path)
    alignment_inspect.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_inspect)


def run_alignment_summary_command(args: Any) -> int | None:
    if args.alignment_command == "alphabet":
        records = load_fasta_alignment(args.alignment)
        alphabet = infer_alignment_alphabet(records)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={"alphabet": alphabet, "sequence_count": len(records)},
                data={
                    "alignment_path": args.alignment,
                    "inferred_alphabet": alphabet,
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "profiles":
        report = list_alignment_filter_profiles()
        outputs = _finalize_outputs(args, command="alignment", inputs=[])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[],
                outputs=outputs,
                metrics={"profile_count": len(report)},
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "gc":
        report = summarise_fasta(args.alignment)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "alphabet": report.inferred_alphabet,
                    "gc_sequence_count": len(report.per_sequence_gc_content),
                },
                data={
                    "alignment_path": report.path,
                    "inferred_alphabet": report.inferred_alphabet,
                    "per_sequence_gc_content": report.per_sequence_gc_content,
                    "whole_alignment_gc_content": report.whole_alignment_gc_content,
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "inspect":
        report = summarise_fasta(args.alignment)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "sequence_count": report.sequence_count,
                    "alignment_length": report.alignment_length,
                    "variable_site_count": report.variable_site_count,
                    "parsimony_informative_site_count": report.parsimony_informative_site_count,
                    "alphabet": report.inferred_alphabet,
                    "invalid_character_count": len(report.invalid_characters),
                    "composition_outlier_count": len(report.composition_outliers),
                    "duplicate_group_count": len(report.duplicate_sequence_groups),
                    "near_duplicate_count": len(report.near_duplicate_pairs),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    return None
