from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.engines import (
    export_workflow_result_bundle,
    inspect_external_engine_preflight,
    replay_workflow_manifest,
    require_preflight_workflow,
    run_phylo_workflow_config,
    validate_workflow_result_bundle,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError
from bijux_phylogenetics.runtime.results import build_command_result


def run_phylo_command(args: Any) -> int:
    executables = {
        "mafft": getattr(args, "mafft_executable", None),
        "trimal": getattr(args, "trimal_executable", None),
        "iqtree": getattr(args, "iqtree_executable", None),
        "fasttree": getattr(args, "fasttree_executable", None),
        "mrbayes": getattr(args, "mrbayes_executable", None),
        "beast": getattr(args, "beast_executable", None),
    }
    if args.phylo_command == "run":
        report = run_phylo_workflow_config(args.config_path)
        outputs = _finalize_outputs(
            args,
            command="phylo",
            inputs=[args.config_path],
            outputs=[
                report.fasta_to_tree_report.manifest_path,
                report.bundle_report.bundle_root,
                report.bundle_report.bundle_manifest_path,
                report.bundle_report.report_path,
            ],
        )
        _print_result(
            build_command_result(
                command="phylo",
                inputs=[args.config_path],
                outputs=outputs,
                warnings=report.warnings + report.notes,
                metrics={
                    "workflow": report.workflow,
                    "selected_workflow_status": (
                        report.selected_workflow_status.readiness_status
                    ),
                    "metadata_present": report.workflow_config.metadata_path is not None,
                    "traits_present": report.workflow_config.traits_path is not None,
                    "alignment_mode": report.workflow_config.alignment_mode,
                    "trimming_mode": report.workflow_config.trimming_mode,
                    "bootstrap_replicates": report.workflow_config.bootstrap_replicates,
                    "iqtree_seed": report.workflow_config.iqtree_seed,
                    "iqtree_threads": report.workflow_config.iqtree_threads,
                    "timeout_seconds": report.workflow_config.timeout_seconds,
                    "bundle_file_count": report.bundle_report.file_count,
                    "bundle_validation_passed": report.bundle_validation.valid,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.phylo_command == "preflight":
        report = inspect_external_engine_preflight(
            executables=executables,
            selected_workflow=args.workflow,
        )
        selected_workflow_status = None
        if args.workflow is not None:
            selected_workflow_status = require_preflight_workflow(
                report, workflow_id=args.workflow
            ).readiness_status
        inputs = [] if args.workflow is None else [args.workflow]
        outputs = _finalize_outputs(args, command="phylo", inputs=inputs)
        _print_result(
            build_command_result(
                command="phylo",
                inputs=inputs,
                outputs=outputs,
                metrics={
                    "engine_count": len(report.engines),
                    "available_engine_count": sum(
                        1 for engine in report.engines if engine.available
                    ),
                    "workflow_count": len(report.workflows),
                    "runnable_workflow_count": sum(
                        1 for workflow in report.workflows if workflow.runnable
                    ),
                    "selected_workflow": args.workflow,
                    "selected_workflow_status": selected_workflow_status,
                    "overall_status": report.overall_status,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.phylo_command == "bundle":
        report = export_workflow_result_bundle(
            args.manifest_path,
            bundle_root=args.out_dir,
        )
        validation = validate_workflow_result_bundle(report.bundle_root)
        outputs = _finalize_outputs(
            args,
            command="phylo",
            inputs=[args.manifest_path],
            outputs=[
                report.bundle_root,
                report.bundle_manifest_path,
                report.report_path,
            ],
        )
        _print_result(
            build_command_result(
                command="phylo",
                inputs=[args.manifest_path],
                outputs=outputs,
                warnings=report.notes,
                metrics={
                    "workflow": report.workflow,
                    "file_count": report.file_count,
                    "copied_input_count": report.copied_input_count,
                    "copied_output_count": report.copied_output_count,
                    "copied_step_manifest_count": report.copied_step_manifest_count,
                    "copied_step_output_count": report.copied_step_output_count,
                    "copied_report_count": report.copied_report_count,
                    "missing_input_count": len(report.missing_input_paths),
                    "validation_passed": validation.valid,
                },
                data={"bundle": report, "validation": validation},
            ),
            json_output=args.json,
        )
        return 0
    if args.phylo_command == "validate-bundle":
        report = validate_workflow_result_bundle(args.bundle_root)
        if not report.valid:
            raise EngineWorkflowError(
                "workflow result bundle validation failed",
                code="workflow_bundle_validation_failed",
                details={
                    "issue_count": len(report.issues),
                    "bundle_root": str(args.bundle_root),
                },
            )
        outputs = _finalize_outputs(
            args,
            command="phylo",
            inputs=[args.bundle_root],
        )
        _print_result(
            build_command_result(
                command="phylo",
                inputs=[args.bundle_root],
                outputs=outputs,
                metrics={
                    "workflow": report.workflow,
                    "file_count": report.file_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    report = replay_workflow_manifest(
        args.manifest_path,
        out_dir=args.out_dir,
        executables=executables,
    )
    outputs = _finalize_outputs(
        args,
        command="phylo",
        inputs=[args.manifest_path],
        outputs=[report.replay_manifest_path],
    )
    _print_result(
        build_command_result(
            command="phylo",
            inputs=[args.manifest_path],
            outputs=outputs,
            metrics={
                "workflow": report.workflow,
                "input_drift_count": len(report.input_drift),
                "changed_input_count": sum(
                    1 for drift in report.input_drift if not drift.matched
                ),
                "engine_version_drift_count": sum(
                    1 for drift in report.engine_version_drift if not drift.matched
                ),
                "comparison_count": len(report.comparisons),
                "outputs_equivalent": report.outputs_equivalent,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0
