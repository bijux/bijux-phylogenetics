from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports.service import (
    render_level_one_release_gate_report,
    render_production_scale_readiness_report,
    render_release_truth_report,
    render_workflow_validation_report,
)
from bijux_phylogenetics.runtime.results import build_command_result


def _build_production_scale_alignment_classes(
    args: Any,
) -> list[tuple[str, int, int]] | None:
    if args.sequence_counts is None and args.alignment_lengths is None:
        return None
    sequence_counts = args.sequence_counts or []
    alignment_lengths = args.alignment_lengths or []
    if len(sequence_counts) != len(alignment_lengths):
        raise ValueError(
            "report production-scale-readiness requires the same number of --sequence-count and --alignment-length values"
        )
    return [
        (
            f"sequences-{sequence_count}-sites-{alignment_length}",
            sequence_count,
            alignment_length,
        )
        for sequence_count, alignment_length in zip(
            sequence_counts,
            alignment_lengths,
            strict=True,
        )
    ]


def _build_production_scale_tree_set_classes(
    args: Any,
) -> list[tuple[str, int, int]] | None:
    if args.posterior_tree_counts is None and args.tree_set_tip_counts is None:
        return None
    posterior_tree_counts = args.posterior_tree_counts or []
    tree_set_tip_counts = args.tree_set_tip_counts or []
    if len(posterior_tree_counts) != len(tree_set_tip_counts):
        raise ValueError(
            "report production-scale-readiness requires the same number of --posterior-tree-count and --tree-set-tip-count values"
        )
    return [
        (f"trees-{tree_count}-taxa-{tip_count}", tree_count, tip_count)
        for tree_count, tip_count in zip(
            posterior_tree_counts,
            tree_set_tip_counts,
            strict=True,
        )
    ]


def add_governance_report_commands(report_subparsers: Any) -> None:
    report_workflow_validation = report_subparsers.add_parser(
        "workflow-validation",
        help="Render the Level 1 workflow validation fixture report.",
    )
    report_workflow_validation.add_argument("--fixtures-root", type=Path)
    report_workflow_validation.add_argument("--out", required=True, type=Path)
    report_workflow_validation.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_workflow_validation)

    report_release_gate = report_subparsers.add_parser(
        "release-gate",
        help="Render the Level 1 release gate for the checked-in workflow fixtures.",
    )
    report_release_gate.add_argument("--fixtures-root", type=Path)
    report_release_gate.add_argument("--out", required=True, type=Path)
    report_release_gate.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_release_gate)

    report_release_truth = report_subparsers.add_parser(
        "release-truth",
        help="Render one machine-produced release truth report from pytest and workflow evidence.",
    )
    report_release_truth.add_argument(
        "--test-report",
        type=Path,
        action="append",
        required=True,
        help="Path to one pytest JUnit XML report for the full test surface. Repeat to aggregate multiple sessions.",
    )
    report_release_truth.add_argument(
        "--real-engine-test-report",
        type=Path,
        action="append",
        required=True,
        help="Path to one pytest JUnit XML report for real-engine tests. Repeat to aggregate multiple sessions.",
    )
    report_release_truth.add_argument("--fixtures-root", type=Path)
    report_release_truth.add_argument(
        "--stress-tier",
        choices=("small", "heavy"),
        default="small",
        help="Governed stress tier to benchmark during release truth generation.",
    )
    report_release_truth.add_argument(
        "--parity-extended",
        action="store_true",
        help="Include the governed extended reference-parity suite in the release truth report.",
    )
    report_release_truth.add_argument("--out", required=True, type=Path)
    report_release_truth.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_release_truth)

    report_production_scale_readiness = report_subparsers.add_parser(
        "production-scale-readiness",
        help="Render one reviewer-facing production-scale readiness report from governed benchmark evidence.",
    )
    report_production_scale_readiness.add_argument("--replicates", type=int, default=1)
    report_production_scale_readiness.add_argument(
        "--tree-tip-count",
        action="append",
        dest="tree_tip_counts",
        type=int,
        help="Add one large-tree taxon count. Repeat to override the governed tree-size classes.",
    )
    report_production_scale_readiness.add_argument(
        "--sequence-count",
        action="append",
        dest="sequence_counts",
        type=int,
        help="Add one sequence count for the large-alignment classes. Repeat alongside --alignment-length.",
    )
    report_production_scale_readiness.add_argument(
        "--alignment-length",
        action="append",
        dest="alignment_lengths",
        type=int,
        help="Add one aligned-site count for the large-alignment classes. Repeat alongside --sequence-count.",
    )
    report_production_scale_readiness.add_argument(
        "--posterior-tree-count",
        action="append",
        dest="posterior_tree_counts",
        type=int,
        help="Add one posterior tree count for the tree-set classes. Repeat alongside --tree-set-tip-count.",
    )
    report_production_scale_readiness.add_argument(
        "--tree-set-tip-count",
        action="append",
        dest="tree_set_tip_counts",
        type=int,
        help="Add one taxon count for the tree-set classes. Repeat alongside --posterior-tree-count.",
    )
    report_production_scale_readiness.add_argument(
        "--stress-tier",
        action="append",
        dest="stress_tiers",
        choices=("small", "heavy"),
        help="Include one governed stress tier. Repeat to aggregate multiple tiers.",
    )
    report_production_scale_readiness.add_argument("--out", required=True, type=Path)
    report_production_scale_readiness.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_production_scale_readiness)


def run_governance_report_command(args: Any) -> int | None:
    if args.report_command == "workflow-validation":
        result = render_workflow_validation_report(
            out_path=args.out,
            fixtures_root=args.fixtures_root,
        )
        inputs = [] if args.fixtures_root is None else [args.fixtures_root]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    metrics={
                        "total_fixture_count": result.validation.total_fixture_count,
                        "passed_fixture_count": result.validation.passed_fixture_count,
                        "workflow_count": len(result.validation.workflows),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "release-gate":
        result = render_level_one_release_gate_report(
            out_path=args.out,
            fixtures_root=args.fixtures_root,
        )
        inputs = [] if args.fixtures_root is None else [args.fixtures_root]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.release_gate.dataset_warnings,
                    metrics={
                        "decision": result.release_gate.gate.decision,
                        "retained_taxa": len(result.release_gate.gate.retained_taxa),
                        "excluded_taxa": len(result.release_gate.gate.excluded_taxa),
                        "blocked_analysis_count": len(
                            result.release_gate.gate.blocked_analyses
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "release-truth":
        result = render_release_truth_report(
            out_path=args.out,
            test_report_paths=args.test_report,
            real_engine_test_report_paths=args.real_engine_test_report,
            fixtures_root=args.fixtures_root,
            include_extended_parity=args.parity_extended,
            stress_tier=args.stress_tier,
        )
        inputs = [
            *args.test_report,
            *args.real_engine_test_report,
            *([args.fixtures_root] if args.fixtures_root is not None else []),
        ]
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=inputs,
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=inputs,
                    outputs=outputs,
                    warnings=result.release_truth.known_limitations,
                    metrics={
                        "total_tests": result.release_truth.total_tests.total_tests,
                        "total_tests_passed": result.release_truth.total_tests.passed_tests,
                        "total_tests_failed": result.release_truth.total_tests.failed_tests,
                        "total_tests_skipped": result.release_truth.total_tests.skipped_tests,
                        "real_engine_tests": result.release_truth.real_engine_tests.total_tests,
                        "real_engine_tests_passed": result.release_truth.real_engine_tests.passed_tests,
                        "real_engine_tests_failed": result.release_truth.real_engine_tests.failed_tests,
                        "real_engine_tests_skipped": result.release_truth.real_engine_tests.skipped_tests,
                        "supported_workflow_count": len(
                            result.release_truth.supported_workflows
                        ),
                        "experimental_workflow_count": len(
                            result.release_truth.experimental_workflows
                        ),
                        "flagship_dataset_count": len(
                            result.release_truth.flagship_datasets
                        ),
                        "reference_parity_case_count": result.release_truth.reference_parity.case_count,
                        "stress_workload_count": len(
                            result.release_truth.stress_suite.observations
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "production-scale-readiness":
        result = render_production_scale_readiness_report(
            out_path=args.out,
            replicates=args.replicates,
            tree_tip_counts=args.tree_tip_counts,
            alignment_size_classes=_build_production_scale_alignment_classes(args),
            tree_set_size_classes=_build_production_scale_tree_set_classes(args),
            stress_tiers=args.stress_tiers,
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[],
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[],
                    outputs=outputs,
                    warnings=result.production_scale_readiness.limitations,
                    metrics=result.machine_manifest["metrics"],
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    return None
