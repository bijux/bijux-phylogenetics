from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.fasta.quality import (
    assess_alignment_low_information,
    build_alignment_forensic_report,
    build_alignment_quality_report,
    build_ambiguous_alignment_column_report,
    build_sequence_quality_ranking,
    detect_over_aligned_regions,
    detect_under_aligned_regions,
    summarize_alignment_readiness,
    summarize_alignment_windows,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_alignment_quality_commands(alignment_subparsers: Any) -> None:
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


def run_alignment_quality_command(args: Any) -> int | None:
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
                    "sequence_length_outlier_count": len(
                        report.sequence_length_outliers
                    ),
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
