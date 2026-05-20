from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.fasta import (
    detect_invalid_alignment_characters,
)
from bijux_phylogenetics.io.fasta.cleaning import (
    detect_composition_outlier_sequences,
    detect_identical_duplicate_sequences,
    detect_near_duplicate_sequences,
)
from bijux_phylogenetics.io.fasta.matrix import (
    compute_alignment_base_frequency_report,
    compute_alignment_segregating_site_report,
    write_alignment_base_frequency_table,
    write_alignment_segregating_site_table,
)
from bijux_phylogenetics.io.fasta.quality import (
    assess_alignment_low_information,
    build_alignment_forensic_report,
    build_alignment_quality_report,
    build_ambiguous_alignment_column_report,
    build_duplicate_sequence_policy_report,
    build_sequence_quality_ranking,
    detect_over_aligned_regions,
    detect_under_aligned_regions,
    summarize_alignment_readiness,
    summarize_alignment_windows,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.runtime.results import build_command_result

from .input_validation import (
    add_alignment_input_validation_commands,
    run_alignment_input_validation_command,
)
from .summary import (
    add_alignment_summary_commands,
    run_alignment_summary_command,
)


def add_alignment_review_commands(alignment_subparsers: Any) -> None:
    add_alignment_summary_commands(alignment_subparsers)
    add_alignment_input_validation_commands(alignment_subparsers)

    alignment_quality = alignment_subparsers.add_parser(
        "quality",
        help="Generate a higher-level alignment quality report.",
    )
    alignment_quality.add_argument("alignment", type=Path)
    alignment_quality.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_quality)

    alignment_windows = alignment_subparsers.add_parser(
        "windows",
        help="Summarize sliding-window alignment quality and suspicious regions.",
    )
    alignment_windows.add_argument("alignment", type=Path)
    alignment_windows.add_argument("--window-size", type=int, default=30)
    alignment_windows.add_argument("--step-size", type=int, default=10)
    alignment_windows.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_windows)

    alignment_readiness = alignment_subparsers.add_parser(
        "readiness",
        help="Classify whether an alignment is ready for distance, ML, Bayesian, coding, or protein workflows.",
    )
    alignment_readiness.add_argument("alignment", type=Path)
    alignment_readiness.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_readiness)

    alignment_forensic = alignment_subparsers.add_parser(
        "forensic",
        help="Build a reviewer-facing alignment forensic report.",
    )
    alignment_forensic.add_argument("alignment", type=Path)
    alignment_forensic.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_forensic)

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

    alignment_duplicates = alignment_subparsers.add_parser(
        "duplicates",
        help="Report identical and near-duplicate aligned sequences.",
    )
    alignment_duplicates.add_argument("alignment", type=Path)
    alignment_duplicates.add_argument("--identity-threshold", type=float, default=0.95)
    alignment_duplicates.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_duplicates)

    alignment_duplicate_policy = alignment_subparsers.add_parser(
        "duplicate-policy",
        help="Recommend how exact and near-duplicate sequences should be handled before inference.",
    )
    alignment_duplicate_policy.add_argument("alignment", type=Path)
    alignment_duplicate_policy.add_argument(
        "--identity-threshold",
        type=float,
        default=0.99,
    )
    alignment_duplicate_policy.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_duplicate_policy)

    alignment_outliers = alignment_subparsers.add_parser(
        "outliers",
        help="Report composition outlier sequences from an alignment.",
    )
    alignment_outliers.add_argument("alignment", type=Path)
    alignment_outliers.add_argument("--deviation-threshold", type=float, default=0.25)
    alignment_outliers.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_outliers)

    alignment_low_information = alignment_subparsers.add_parser(
        "low-information",
        help="Report whether an alignment has enough informative sites for defensible inference.",
    )
    alignment_low_information.add_argument("alignment", type=Path)
    alignment_low_information.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_low_information)

    alignment_ambiguous_columns = alignment_subparsers.add_parser(
        "ambiguous-columns",
        help="List columns dominated by ambiguity, missing data, or gaps.",
    )
    alignment_ambiguous_columns.add_argument("alignment", type=Path)
    alignment_ambiguous_columns.add_argument("--threshold", type=float, default=0.5)
    alignment_ambiguous_columns.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_ambiguous_columns)

    alignment_sequence_ranking = alignment_subparsers.add_parser(
        "sequence-ranking",
        help="Rank sequences by missingness, ambiguity, gap burden, composition, and duplicate status.",
    )
    alignment_sequence_ranking.add_argument("alignment", type=Path)
    alignment_sequence_ranking.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(alignment_sequence_ranking)


def run_alignment_review_command(args: Any) -> int | None:
    summary_exit_code = run_alignment_summary_command(args)
    if summary_exit_code is not None:
        return summary_exit_code

    input_validation_exit_code = run_alignment_input_validation_command(args)
    if input_validation_exit_code is not None:
        return input_validation_exit_code
    if args.alignment_command == "quality":
        report = build_alignment_quality_report(args.alignment)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "quality_score": report.quality_score,
                    "invariant_site_count": report.invariant_site_count,
                    "parsimony_informative_site_count": report.parsimony_informative_site_count,
                    "suspicious_alignment": report.suspicious_alignment,
                    "suspicious_reason_count": len(report.suspicious_reasons),
                    "concentrated_column_count": report.missing_data_concentration.concentrated_column_count,
                    "invalid_character_count": len(report.invalid_characters),
                    "composition_outlier_count": len(report.composition_outliers),
                    "sequence_length_outlier_count": len(report.sequence_length_outliers),
                    "duplicate_group_count": len(report.duplicate_sequence_groups),
                    "near_duplicate_count": len(report.near_duplicate_pairs),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "windows":
        windows = summarize_alignment_windows(
            args.alignment,
            window_size=args.window_size,
            step_size=args.step_size,
        )
        over_aligned = detect_over_aligned_regions(
            args.alignment,
            window_size=args.window_size,
            step_size=args.step_size,
        )
        under_aligned = detect_under_aligned_regions(
            args.alignment,
            window_size=args.window_size,
            step_size=args.step_size,
        )
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=[region.note for region in over_aligned + under_aligned],
                metrics={
                    "window_count": len(windows),
                    "over_aligned_region_count": len(over_aligned),
                    "under_aligned_region_count": len(under_aligned),
                },
                data={
                    "windows": windows,
                    "over_aligned_regions": over_aligned,
                    "under_aligned_regions": under_aligned,
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "readiness":
        report = summarize_alignment_readiness(args.alignment)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "sequence_count": report.sequence_count,
                    "alignment_length": report.alignment_length,
                    "ready_method_count": sum(
                        1 for method in report.methods if method.ready
                    ),
                    "blocked_method_count": sum(
                        1 for method in report.methods if not method.ready
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "forensic":
        report = build_alignment_forensic_report(args.alignment)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "quality_score": report.quality.quality_score,
                    "safe_for_distance_analysis": report.safe_for_distance_analysis,
                    "safe_for_maximum_likelihood": report.safe_for_maximum_likelihood,
                    "safe_for_bayesian_inference": report.safe_for_bayesian_inference,
                    "safe_for_coding_analysis": report.safe_for_coding_analysis,
                    "safe_for_publication": report.safe_for_publication,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
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
    if args.alignment_command == "duplicates":
        duplicates = detect_identical_duplicate_sequences(args.alignment)
        near_duplicates = detect_near_duplicate_sequences(
            args.alignment,
            identity_threshold=args.identity_threshold,
        )
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "duplicate_group_count": len(duplicates),
                    "near_duplicate_count": len(near_duplicates),
                    "identity_threshold": args.identity_threshold,
                },
                data={
                    "duplicate_sequence_groups": duplicates,
                    "near_duplicate_pairs": near_duplicates,
                },
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "duplicate-policy":
        report = build_duplicate_sequence_policy_report(
            args.alignment,
            near_duplicate_threshold=args.identity_threshold,
        )
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "exact_duplicate_group_count": len(report.exact_duplicate_groups),
                    "near_duplicate_pair_count": len(report.near_duplicate_pairs),
                    "policy_action_count": len(report.policy_actions),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "outliers":
        report = detect_composition_outlier_sequences(
            args.alignment,
            deviation_threshold=args.deviation_threshold,
        )
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "composition_outlier_count": len(report),
                    "deviation_threshold": args.deviation_threshold,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "low-information":
        report = assess_alignment_low_information(args.alignment)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.reasons,
                metrics={
                    "low_information": report.low_information,
                    "parsimony_informative_site_count": report.parsimony_informative_site_count,
                    "alignment_length": report.alignment_length,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "ambiguous-columns":
        report = build_ambiguous_alignment_column_report(
            args.alignment,
            threshold=args.threshold,
        )
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "ambiguous_column_count": len(report.rows),
                    "threshold": report.threshold,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.alignment_command == "sequence-ranking":
        report = build_sequence_quality_ranking(args.alignment)
        outputs = _finalize_outputs(args, command="alignment", inputs=[args.alignment])
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "sequence_count": len(report.rows),
                    "lowest_score": None if not report.rows else report.rows[0].score,
                    "highest_score": None if not report.rows else report.rows[-1].score,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    return None
