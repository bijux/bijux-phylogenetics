from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    render_calibration_audit_report,
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
)
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    method_tier_metrics,
    method_tier_warnings,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_beast_calibration_review_commands(adapter_subparsers: Any) -> None:
    adapter_beast_calibrations = adapter_subparsers.add_parser(
        "beast-calibrations",
        help="Validate a fossil calibration table against a tree.",
    )
    adapter_beast_calibrations.add_argument("tree_path", type=Path)
    adapter_beast_calibrations.add_argument("calibration_path", type=Path)
    adapter_beast_calibrations.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(adapter_beast_calibrations)

    adapter_beast_tip_dates = adapter_subparsers.add_parser(
        "beast-tip-dates",
        help="Validate tip-dating metadata against a tree and optional alignment.",
    )
    adapter_beast_tip_dates.add_argument("tree_path", type=Path)
    adapter_beast_tip_dates.add_argument("tip_dates_path", type=Path)
    adapter_beast_tip_dates.add_argument("--alignment", type=Path)
    adapter_beast_tip_dates.add_argument("--date-column", default="date")
    adapter_beast_tip_dates.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(adapter_beast_tip_dates)

    adapter_beast_calibration_report = adapter_subparsers.add_parser(
        "beast-calibration-report",
        help="Render an HTML calibration audit report.",
    )
    adapter_beast_calibration_report.add_argument("tree_path", type=Path)
    adapter_beast_calibration_report.add_argument("calibration_path", type=Path)
    adapter_beast_calibration_report.add_argument("--out", required=True, type=Path)
    adapter_beast_calibration_report.add_argument("--tip-dates", type=Path)
    adapter_beast_calibration_report.add_argument("--alignment", type=Path)
    adapter_beast_calibration_report.add_argument("--date-column", default="date")
    adapter_beast_calibration_report.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(adapter_beast_calibration_report)


def run_beast_calibration_review_command(args: Any) -> int | None:
    if args.adapter_command == "beast-calibrations":
        report = validate_fossil_calibration_table(
            args.tree_path, args.calibration_path
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.tree_path, args.calibration_path],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.tree_path, args.calibration_path],
                outputs=outputs,
                metrics={
                    "calibration_count": report.calibration_count,
                    "invalid_calibration_count": report.invalid_calibration_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "beast-tip-dates":
        report = validate_tip_dating_metadata(
            args.tree_path,
            args.tip_dates_path,
            alignment_path=args.alignment,
            date_column=args.date_column,
        )
        inputs = [
            args.tree_path,
            args.tip_dates_path,
            *([args.alignment] if args.alignment is not None else []),
        ]
        outputs = _finalize_outputs(args, command="adapter", inputs=inputs)
        _print_result(
            build_command_result(
                command="adapter",
                inputs=inputs,
                outputs=outputs,
                metrics={
                    "valid_tip_count": report.valid_tip_count,
                    "invalid_tip_count": report.invalid_tip_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command != "beast-calibration-report":
        return None

    report = render_calibration_audit_report(
        tree_path=args.tree_path,
        calibration_path=args.calibration_path,
        out_path=args.out,
        tip_dates_path=args.tip_dates,
        alignment_path=args.alignment,
        date_column=args.date_column,
    )
    inputs = [
        args.tree_path,
        args.calibration_path,
        *([args.tip_dates] if args.tip_dates is not None else []),
        *([args.alignment] if args.alignment is not None else []),
    ]
    outputs = _finalize_outputs(
        args,
        command="adapter",
        inputs=inputs,
        outputs=[report.output_path],
    )
    _print_result(
        build_command_result(
            command="adapter",
            inputs=inputs,
            outputs=outputs,
            warnings=method_tier_warnings(report.method_tier),
            metrics={
                "invalid_calibration_count": report.invalid_calibration_count,
                "warning_count": report.warning_count
                + len(method_tier_warnings(report.method_tier)),
                **method_tier_metrics(report.method_tier),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
