from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.evidence.provenance.method_tiers import method_tier_metrics
from bijux_phylogenetics.reports import (
    write_publication_package_comparison_report,
    write_publication_package_revalidation_report,
    write_reviewer_audit_checklist_from_manifest,
)
from bijux_phylogenetics.reports.publication.tree import build_tree_report_package
from bijux_phylogenetics.reports.service import render_tree_report
from bijux_phylogenetics.runtime.results import build_command_result


def add_publication_report_commands(report_subparsers: Any) -> None:
    report_tree = report_subparsers.add_parser(
        "tree", help="Render a deterministic single-tree HTML report."
    )
    report_tree.add_argument("tree", type=Path)
    report_tree.add_argument("--out", required=True, type=Path)
    report_tree.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_tree)

    report_tree_package = report_subparsers.add_parser(
        "tree-package",
        help="Build a full tree report package with figure and TSV ledgers.",
    )
    report_tree_package.add_argument("tree", type=Path)
    report_tree_package.add_argument("--out-dir", required=True, type=Path)
    report_tree_package.add_argument(
        "--json", action="store_true", help="Emit the package build result as JSON."
    )
    _add_manifest_argument(report_tree_package)

    report_reviewer_audit_checklist = report_subparsers.add_parser(
        "reviewer-audit-checklist",
        help="Write a reviewer-facing audit checklist from a supported package manifest.",
    )
    report_reviewer_audit_checklist.add_argument("manifest", type=Path)
    report_reviewer_audit_checklist.add_argument("--out", required=True, type=Path)
    report_reviewer_audit_checklist.add_argument(
        "--json", action="store_true", help="Emit the checklist result as JSON."
    )
    _add_manifest_argument(report_reviewer_audit_checklist)

    report_package_revalidation = report_subparsers.add_parser(
        "package-revalidation",
        help="Revalidate a stored publication package from its manifest and checksum inventory.",
    )
    report_package_revalidation.add_argument("manifest", type=Path)
    report_package_revalidation.add_argument("--out-dir", required=True, type=Path)
    report_package_revalidation.add_argument(
        "--json",
        action="store_true",
        help="Emit the revalidation result as JSON.",
    )
    _add_manifest_argument(report_package_revalidation)

    report_package_comparison = report_subparsers.add_parser(
        "package-comparison",
        help="Compare two stored publication package versions for the same governed study.",
    )
    report_package_comparison.add_argument("left_manifest", type=Path)
    report_package_comparison.add_argument("right_manifest", type=Path)
    report_package_comparison.add_argument("--out-dir", required=True, type=Path)
    report_package_comparison.add_argument(
        "--json",
        action="store_true",
        help="Emit the comparison result as JSON.",
    )
    _add_manifest_argument(report_package_comparison)


def run_publication_report_command(args: Any) -> int | None:
    if args.report_command == "tree":
        result = render_tree_report(tree_path=args.tree, out_path=args.out)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.tree],
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.tree],
                    outputs=outputs,
                    warnings=result.validation.warnings + result.inspection.warnings,
                    metrics={"tip_count": result.inspection.tip_count},
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "tree-package":
        result = build_tree_report_package(args.tree, out_dir=args.out_dir)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.tree],
            outputs=[
                result.report_path,
                result.figure_path,
                result.methods_summary_path,
                result.reviewer_audit_checklist_path,
                result.support_table_path,
                result.clade_table_path,
                result.branch_stats_path,
                result.manifest_path,
            ],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.tree],
                    outputs=outputs,
                    warnings=result.validation.warnings + result.inspection.warnings,
                    metrics={
                        "tip_count": result.inspection.tip_count,
                        "supported_branch_count": sum(
                            1 for row in result.support_rows if row.support is not None
                        ),
                        "rendered_support_count": result.figure.rendered_support_count,
                        "long_outlier_count": result.branch_stats.long_outlier_count,
                        "methods_warning_count": result.methods_summary.warning_count,
                        "reviewer_audit_item_count": len(
                            result.reviewer_audit_checklist.items
                        ),
                        **method_tier_metrics(result.method_tier),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_dir)
        return 0

    if args.report_command == "reviewer-audit-checklist":
        result = write_reviewer_audit_checklist_from_manifest(args.out, args.manifest)
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.manifest],
            outputs=[result.output_path],
        )
        if args.json:
            blocked_item_count = sum(
                1 for item in result.checklist.items if item.status == "blocked"
            )
            risk_item_count = sum(
                1 for item in result.checklist.items if item.status == "risk"
            )
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.manifest],
                    outputs=outputs,
                    warnings=[
                        item.summary
                        for item in result.checklist.items
                        if item.status != "pass"
                    ],
                    metrics={
                        "report_kind": result.checklist.report_kind,
                        "item_count": len(result.checklist.items),
                        "blocked_item_count": blocked_item_count,
                        "risk_item_count": risk_item_count,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "package-revalidation":
        result = write_publication_package_revalidation_report(
            args.out_dir,
            args.manifest,
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.manifest],
            outputs=[
                result.artifact_table_path,
                result.check_table_path,
                result.summary_path,
                result.report_path,
            ],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.manifest],
                    outputs=outputs,
                    warnings=[
                        row.summary for row in result.check_rows if row.status != "pass"
                    ],
                    metrics={
                        "report_kind": result.report_kind,
                        "artifact_row_count": len(result.artifact_rows),
                        "check_row_count": len(result.check_rows),
                        "matched_artifact_count": result.matched_artifact_count,
                        "missing_artifact_count": result.missing_artifact_count,
                        "checksum_mismatch_count": result.checksum_mismatch_count,
                        "size_mismatch_count": result.size_mismatch_count,
                        "unexpected_file_count": result.unexpected_file_count,
                        "blocked_check_count": result.blocked_check_count,
                        "risk_check_count": result.risk_check_count,
                        "all_original_artifacts_match": (
                            result.all_original_artifacts_match
                        ),
                        "overall_revalidation_status": (
                            result.overall_revalidation_status
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    if args.report_command == "package-comparison":
        result = write_publication_package_comparison_report(
            args.out_dir,
            args.left_manifest,
            args.right_manifest,
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.left_manifest, args.right_manifest],
            outputs=[
                result.artifact_table_path,
                result.check_table_path,
                result.summary_path,
                result.report_path,
            ],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.left_manifest, args.right_manifest],
                    outputs=outputs,
                    warnings=[
                        row.summary for row in result.check_rows if row.status != "pass"
                    ],
                    metrics={
                        "report_kind": result.report_kind,
                        "dataset_id": result.dataset_id,
                        "artifact_row_count": len(result.artifact_rows),
                        "check_row_count": len(result.check_rows),
                        "same_artifact_count": result.same_artifact_count,
                        "changed_artifact_count": result.changed_artifact_count,
                        "left_only_artifact_count": result.left_only_artifact_count,
                        "right_only_artifact_count": result.right_only_artifact_count,
                        "config_difference_count": result.config_difference_count,
                        "sequence_left_only_count": result.sequence_left_only_count,
                        "sequence_right_only_count": result.sequence_right_only_count,
                        "accession_left_only_count": result.accession_left_only_count,
                        "accession_right_only_count": result.accession_right_only_count,
                        "alignment_difference_count": result.alignment_difference_count,
                        "figure_or_report_difference_count": (
                            result.figure_or_report_difference_count
                        ),
                        "scientific_finding_difference_count": (
                            result.scientific_finding_difference_count
                        ),
                        "overall_comparison_status": (result.overall_comparison_status),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_root)
        return 0

    return None
