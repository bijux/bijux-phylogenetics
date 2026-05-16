from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.provenance.method_tiers import method_tier_metrics
from bijux_phylogenetics.reports.service import (
    render_alignment_report,
    render_dataset_report,
    render_level_one_release_gate_report,
    render_phylo_inputs_report,
    render_release_truth_report,
    render_taxon_report,
    render_tree_report,
    render_workflow_validation_report,
)
from bijux_phylogenetics.reports.tree_package import build_tree_report_package
from bijux_phylogenetics.runtime.results import build_command_result


def add_report_command(subparsers: Any) -> None:
    report = subparsers.add_parser(
        get_command_spec("report").name, help=get_command_spec("report").summary
    )
    report_subparsers = report.add_subparsers(dest="report_command", required=True)

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

    report_alignment = report_subparsers.add_parser(
        "alignment", help="Render an alignment-only HTML diagnostic report."
    )
    report_alignment.add_argument("--alignment", required=True, type=Path)
    report_alignment.add_argument("--out", required=True, type=Path)
    report_alignment.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_alignment)

    report_dataset = report_subparsers.add_parser(
        "dataset", help="Render a tree plus table dataset HTML report."
    )
    report_dataset.add_argument("--tree", required=True, type=Path)
    report_dataset.add_argument("--metadata", required=True, type=Path)
    report_dataset.add_argument("--traits", type=Path)
    report_dataset.add_argument("--alignment", type=Path)
    report_dataset.add_argument("--tip-dates", type=Path)
    report_dataset.add_argument("--calibrations", type=Path)
    report_dataset.add_argument("--out", required=True, type=Path)
    report_dataset.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_dataset)

    report_phylo_inputs = report_subparsers.add_parser(
        "phylo-inputs",
        help="Render a tree plus alignment HTML input report.",
    )
    report_phylo_inputs.add_argument("--tree", required=True, type=Path)
    report_phylo_inputs.add_argument("--alignment", required=True, type=Path)
    report_phylo_inputs.add_argument("--out", required=True, type=Path)
    report_phylo_inputs.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_phylo_inputs)

    report_taxonomy = report_subparsers.add_parser(
        "taxonomy", help="Render a reviewer-facing taxon audit HTML report."
    )
    report_taxonomy.add_argument("--tree", required=True, type=Path)
    report_taxonomy.add_argument("--synonym-table", type=Path)
    report_taxonomy.add_argument("--metadata", type=Path)
    report_taxonomy.add_argument("--traits", type=Path)
    report_taxonomy.add_argument("--alignment", type=Path)
    report_taxonomy.add_argument("--filtered-alignment", type=Path)
    report_taxonomy.add_argument("--inference-tree", type=Path)
    report_taxonomy.add_argument("--reported-taxa", type=Path)
    report_taxonomy.add_argument("--out", required=True, type=Path)
    report_taxonomy.add_argument(
        "--json", action="store_true", help="Emit the report build result as JSON."
    )
    _add_manifest_argument(report_taxonomy)

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


def run_report_command(args: Any) -> int:
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
                        **method_tier_metrics(result.method_tier),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_dir)
        return 0

    if args.report_command == "alignment":
        result = render_alignment_report(
            alignment_path=args.alignment, out_path=args.out
        )
        outputs = _finalize_outputs(
            args,
            command="report",
            inputs=[args.alignment],
            outputs=[result.output_path, result.machine_manifest_path],
        )
        if args.json:
            _print_result(
                build_command_result(
                    command="report",
                    inputs=[args.alignment],
                    outputs=outputs,
                    warnings=result.alignment_forensic.warnings,
                    metrics={
                        "sequence_count": result.alignment.sequence_count,
                        "alignment_length": result.alignment.alignment_length,
                        "quality_score": result.alignment_quality.quality_score,
                        "warning_count": len(result.alignment_forensic.warnings),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "dataset":
        result = render_dataset_report(
            tree_path=args.tree,
            metadata_path=args.metadata,
            traits_path=args.traits,
            alignment_path=args.alignment,
            tip_dates_path=args.tip_dates,
            calibration_path=args.calibrations,
            out_path=args.out,
        )
        inputs = [args.tree, args.metadata]
        if args.traits is not None:
            inputs.append(args.traits)
        if args.alignment is not None:
            inputs.append(args.alignment)
        if args.tip_dates is not None:
            inputs.append(args.tip_dates)
        if args.calibrations is not None:
            inputs.append(args.calibrations)
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
                    warnings=result.validation.warnings + result.inspection.warnings,
                    metrics={
                        "tip_count": result.inspection.tip_count,
                        "linked_taxa": result.metadata_linkage.linked_taxa,
                        "readiness_decision": None
                        if result.dataset_audit is None
                        else result.dataset_audit.readiness_decision,
                        "excluded_taxa": 0
                        if result.dataset_audit is None
                        else len(result.dataset_audit.exclusion_table.rows),
                        "blocked_analysis_count": 0
                        if result.dataset_audit is None
                        else len(result.dataset_audit.blocked_analyses),
                        "risky_analysis_count": 0
                        if result.dataset_audit is None
                        else sum(
                            1
                            for row in result.dataset_audit.analysis_decisions
                            if row.decision == "risky"
                        ),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "phylo-inputs":
        result = render_phylo_inputs_report(
            tree_path=args.tree,
            alignment_path=args.alignment,
            out_path=args.out,
        )
        inputs = [args.tree, args.alignment]
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
                    warnings=result.validation.warnings + result.inspection.warnings,
                    metrics={
                        "tip_count": result.inspection.tip_count,
                        "alignment_length": result.alignment.alignment_length,
                        "linked_taxa": result.alignment_linkage.linked_taxa,
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

    if args.report_command == "taxonomy":
        result = render_taxon_report(
            tree_path=args.tree,
            synonym_table_path=args.synonym_table,
            metadata_path=args.metadata,
            traits_path=args.traits,
            alignment_path=args.alignment,
            filtered_alignment_path=args.filtered_alignment,
            inference_tree_path=args.inference_tree,
            reported_taxa_path=args.reported_taxa,
            out_path=args.out,
        )
        inputs = [args.tree, *([args.synonym_table] if args.synonym_table is not None else [])]
        if args.metadata is not None:
            inputs.append(args.metadata)
        if args.traits is not None:
            inputs.append(args.traits)
        if args.alignment is not None:
            inputs.append(args.alignment)
        if args.filtered_alignment is not None:
            inputs.append(args.filtered_alignment)
        if args.inference_tree is not None:
            inputs.append(args.inference_tree)
        if args.reported_taxa is not None:
            inputs.append(args.reported_taxa)
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
                    warnings=result.taxon_audit.warnings,
                    metrics={
                        "tree_tip_count": result.taxon_audit.tree_tip_count,
                        "status": result.taxon_audit.status,
                        "conflict_count": len(result.taxon_audit.mapping_conflicts.rows),
                        "crosswalk_rows": 0
                        if result.taxon_crosswalk is None
                        else len(result.taxon_crosswalk.rows),
                        "excluded_taxa": 0
                        if result.taxon_exclusions is None
                        else len(result.taxon_exclusions.rows),
                        "loss_stage_count": 0
                        if result.taxon_workflow_loss is None
                        else len(result.taxon_workflow_loss.loss_stage_counts),
                        "unstable_taxa": 0
                        if result.taxon_stability is None
                        else len(result.taxon_stability.unstable_taxa),
                    },
                    data=result,
                ),
                json_output=True,
            )
            return 0
        print(result.output_path)
        return 0

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

    raise NotImplementedError(f"unsupported report command: {args.report_command}")
