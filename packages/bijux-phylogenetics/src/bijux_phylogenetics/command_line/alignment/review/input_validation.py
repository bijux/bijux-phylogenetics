from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.fasta import (
    detect_fasta_sequence_type,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.records import (
    classify_alignment_sequences,
    detect_sequence_length_outliers,
    repair_fasta_input,
    validate_fasta_input,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_alignment_input_validation_commands(alignment_subparsers: Any) -> None:
    alignment_classify = alignment_subparsers.add_parser(
        "classify",
        help="Classify whether a FASTA input is aligned, raw, or shape-ambiguous.",
    )
    alignment_classify.add_argument("alignment", type=Path)
    alignment_classify.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_classify)

    alignment_sequence_type = alignment_subparsers.add_parser(
        "sequence-type",
        help="Classify raw FASTA sequence type compatibility and confidence.",
    )
    alignment_sequence_type.add_argument("alignment", type=Path)
    alignment_sequence_type.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_sequence_type)

    alignment_validate_input = alignment_subparsers.add_parser(
        "validate-input",
        help="Validate raw FASTA input for duplicates, illegal characters, empty sequences, and length outliers.",
    )
    alignment_validate_input.add_argument("alignment", type=Path)
    alignment_validate_input.add_argument(
        "--sequence-type",
        choices=("dna", "rna", "protein"),
    )
    alignment_validate_input.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_validate_input)

    alignment_repair_input = alignment_subparsers.add_parser(
        "repair-input",
        help="Write a repaired FASTA input after explicit identifier normalization or invalid-record removal.",
    )
    alignment_repair_input.add_argument("alignment", type=Path)
    alignment_repair_input.add_argument("--out", required=True, type=Path)
    alignment_repair_input.add_argument(
        "--sequence-type",
        choices=("dna", "rna", "protein"),
    )
    alignment_repair_input.add_argument(
        "--normalize-identifiers",
        action="store_true",
        help="Rewrite FASTA identifiers into engine-safe stable names and resolve collisions.",
    )
    alignment_repair_input.add_argument(
        "--remove-invalid-records",
        action="store_true",
        help="Remove records with empty sequences or unsupported characters.",
    )
    alignment_repair_input.add_argument(
        "--json", action="store_true", help="Emit the repair report as JSON."
    )
    _add_manifest_argument(alignment_repair_input)

    alignment_length_outliers = alignment_subparsers.add_parser(
        "length-outliers",
        help="Report raw sequence length outliers before alignment assumptions are imposed.",
    )
    alignment_length_outliers.add_argument("alignment", type=Path)
    alignment_length_outliers.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_length_outliers)


def run_alignment_input_validation_command(args: Any) -> int | None:
    if args.alignment_command == "classify":
        report = classify_alignment_sequences(args.alignment)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "sequence_count": report.sequence_count,
                    "min_sequence_length": report.min_sequence_length,
                    "max_sequence_length": report.max_sequence_length,
                    "state": report.state,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "sequence-type":
        report = detect_fasta_sequence_type(args.alignment)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "detected_type": report.detected_type,
                    "selected_type": report.selected_type,
                    "confidence": report.confidence,
                    "compatible_type_count": len(report.compatible_types),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "validate-input":
        report = validate_fasta_input(
            args.alignment,
            sequence_type=args.sequence_type,
        )
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "sequence_count": report.summary.sequence_count,
                    "duplicate_identifier_count": len(report.duplicate_identifiers),
                    "illegal_character_count": len(report.illegal_characters),
                    "empty_sequence_count": len(report.empty_sequences),
                    "sequence_length_outlier_count": len(report.length_outliers),
                    "inferred_alphabet": report.summary.inferred_alphabet,
                    "detected_type": report.sequence_type_report.detected_type,
                    "selected_type": report.sequence_type_report.selected_type,
                    "sequence_type_confidence": report.sequence_type_report.confidence,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "repair-input":
        if not args.normalize_identifiers and not args.remove_invalid_records:
            raise ValueError(
                "repair-input requires at least one explicit repair action"
            )
        records, report = repair_fasta_input(
            args.alignment,
            sequence_type=args.sequence_type,
            normalize_identifiers=args.normalize_identifiers,
            remove_invalid_records=args.remove_invalid_records,
        )
        output_path = write_fasta_alignment(args.out, records)
        report.output_path = output_path
        repaired_validation = validate_fasta_input(
            output_path,
            sequence_type=args.sequence_type,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
            outputs=[output_path],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=list(
                    dict.fromkeys(report.warnings + repaired_validation.warnings)
                ),
                metrics={
                    "before_sequence_count": report.before.sequence_count,
                    "after_sequence_count": report.after.sequence_count,
                    "normalized_identifier_count": len(report.normalized_identifiers),
                    "removed_record_count": len(report.removed_records),
                    "remaining_duplicate_identifier_count": len(
                        repaired_validation.duplicate_identifiers
                    ),
                    "remaining_illegal_character_count": len(
                        repaired_validation.illegal_characters
                    ),
                    "remaining_empty_sequence_count": len(
                        repaired_validation.empty_sequences
                    ),
                },
                data={
                    "repair": report,
                    "post_repair_validation": repaired_validation,
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "length-outliers":
        report = detect_sequence_length_outliers(args.alignment)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={"sequence_length_outlier_count": len(report)},
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    return None
