from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_distance_tree_method_argument,
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.distance import (
    assess_distance_method_maturity,
    build_distance_method_report,
    write_distance_reproducibility_bundle,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .shared import (
    add_ambiguity_policy_option,
    add_distance_model_option,
    add_gap_handling_option,
)


def add_distance_reporting_commands(alignment_subparsers: Any) -> None:
    alignment_distance_method_report = alignment_subparsers.add_parser(
        "distance-method-report",
        help="Build a structured distance-method report with support, model, and gap-sensitivity sections.",
    )
    alignment_distance_method_report.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_distance_method_report)
    add_distance_model_option(alignment_distance_method_report)
    add_gap_handling_option(alignment_distance_method_report)
    add_ambiguity_policy_option(alignment_distance_method_report)
    alignment_distance_method_report.add_argument("--replicates", type=int, default=25)
    alignment_distance_method_report.add_argument("--seed", type=int, default=1)
    alignment_distance_method_report.add_argument(
        "--json", action="store_true", help="Emit the structured report as JSON."
    )
    _add_manifest_argument(alignment_distance_method_report)

    alignment_distance_maturity = alignment_subparsers.add_parser(
        "distance-maturity",
        help="Run the distance-method maturity gate for one alignment.",
    )
    alignment_distance_maturity.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_distance_maturity)
    add_distance_model_option(alignment_distance_maturity)
    add_gap_handling_option(alignment_distance_maturity)
    add_ambiguity_policy_option(alignment_distance_maturity)
    alignment_distance_maturity.add_argument("--replicates", type=int, default=25)
    alignment_distance_maturity.add_argument("--seed", type=int, default=1)
    alignment_distance_maturity.add_argument(
        "--json", action="store_true", help="Emit the maturity gate as JSON."
    )
    _add_manifest_argument(alignment_distance_maturity)

    alignment_distance_bundle = alignment_subparsers.add_parser(
        "distance-bundle",
        help="Write a reproducibility bundle for one distance-analysis workflow.",
    )
    alignment_distance_bundle.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_distance_bundle)
    add_distance_model_option(alignment_distance_bundle)
    add_gap_handling_option(alignment_distance_bundle)
    add_ambiguity_policy_option(alignment_distance_bundle)
    alignment_distance_bundle.add_argument("--replicates", type=int, default=100)
    alignment_distance_bundle.add_argument("--seed", type=int, default=1)
    alignment_distance_bundle.add_argument("--out-dir", required=True, type=Path)
    alignment_distance_bundle.add_argument(
        "--json", action="store_true", help="Emit the bundle report as JSON."
    )
    _add_manifest_argument(alignment_distance_bundle)


def run_distance_reporting_command(args: Any) -> int | None:
    if args.alignment_command == "distance-method-report":
        report = build_distance_method_report(
            args.alignment,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
            bootstrap_replicates=args.replicates,
            bootstrap_seed=args.seed,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.maturity_gate.warnings,
                metrics={
                    "method": report.method,
                    "decision": report.maturity_gate.decision,
                    "bootstrap_clade_count": report.bootstrap_summary.clade_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.alignment_command == "distance-maturity":
        report = assess_distance_method_maturity(
            args.alignment,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
            bootstrap_replicates=args.replicates,
            bootstrap_seed=args.seed,
            validate_bundle=True,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "decision": report.decision,
                    "check_count": len(report.checks),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.alignment_command == "distance-bundle":
        report = write_distance_reproducibility_bundle(
            args.out_dir,
            alignment_path=args.alignment,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
            replicates=args.replicates,
            seed=args.seed,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
            outputs=list(report.files),
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "file_count": len(report.files),
                    "replicates": report.replicates,
                    "method": report.method,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None
