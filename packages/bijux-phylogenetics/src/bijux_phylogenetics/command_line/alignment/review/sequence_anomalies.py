from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.io.fasta.cleaning import (
    detect_composition_outlier_sequences,
    detect_identical_duplicate_sequences,
    detect_near_duplicate_sequences,
)
from bijux_phylogenetics.io.fasta.quality import build_duplicate_sequence_policy_report
from bijux_phylogenetics.runtime.results import build_command_result


def add_alignment_sequence_anomaly_commands(alignment_subparsers: Any) -> None:
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


def run_alignment_sequence_anomaly_command(args: Any) -> int | None:
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
    return None
