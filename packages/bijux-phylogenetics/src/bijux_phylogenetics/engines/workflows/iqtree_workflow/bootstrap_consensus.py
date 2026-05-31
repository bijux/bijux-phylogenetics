from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from ...common import (
    build_file_checksums,
    execute_engine_command,
    read_engine_version,
    resolve_engine_executable,
    validate_timeout_seconds,
)
from ...validation.preflight import require_external_engine_surface
from ..models import EngineWorkflowReport
from ..state import (
    _persist_workflow_report,
    _prefix_path,
    _record_output_validation_failure,
    _require_nonempty_text_output,
    _resolve_incomplete_workflow_state,
    _resume_existing_workflow,
    _validate_complete_support_coverage,
    _validate_matching_tree_taxa,
    _validate_support_value_count,
    _validate_tree_output,
)
from .shared import (
    _build_iqtree_model_selection_summary,
    _build_iqtree_summary,
    _existing_iqtree_outputs,
)


def run_bootstrap_consensus_tree(
    bootstrap_trees_path: Path,
    *,
    out_dir: Path,
    prefix: str = "bootstrap-consensus",
    executable: str | Path = "iqtree2",
    minimum_support: float = 0.5,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Construct a consensus tree from bootstrap trees."""
    if not 0.0 <= minimum_support <= 1.0:
        raise ValueError(
            f"minimum_support must be between 0 and 1 inclusive, got {minimum_support}"
        )
    if not bootstrap_trees_path.exists():
        raise FileNotFoundError(bootstrap_trees_path)
    validate_timeout_seconds(timeout_seconds)
    require_external_engine_surface(
        workflow_id="bootstrap-consensus",
        summary="IQ-TREE bootstrap-consensus workflow.",
        required_engines=("iqtree",),
        executables={"iqtree": executable},
        preserve_missing_error=True,
    )
    prefix_path = _prefix_path(out_dir, prefix)
    manifest_path = prefix_path.with_suffix(".manifest.json")
    version = read_engine_version(
        "iqtree",
        executable,
        version_args=("--version",),
        timeout_seconds=timeout_seconds,
    )
    resolved = resolve_engine_executable(executable)
    consensus_tree_path = prefix_path.with_suffix(".contree")
    log_path = prefix_path.with_suffix(".log")
    command = [
        resolved,
        "-t",
        str(bootstrap_trees_path.resolve()),
        "-con",
        "-minsup",
        str(minimum_support),
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[bootstrap_trees_path],
            expected_command=command,
            expected_version=version,
        )
        if resumed is not None:
            return resumed
    incomplete_notes = _resolve_incomplete_workflow_state(
        manifest_path=manifest_path,
        incomplete_run_policy=incomplete_run_policy,
    )
    run = execute_engine_command(
        engine_name="iqtree",
        workflow="bootstrap-consensus",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "consensus_tree": consensus_tree_path,
            "iqtree_log": log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _require_nonempty_text_output(
            log_path,
            engine_name="iqtree",
            workflow="bootstrap-consensus",
            output_name="iqtree_log",
            artifact_kind="iqtree-log",
        )
        _validate_tree_output(
            consensus_tree_path,
            engine_name="iqtree",
            workflow="bootstrap-consensus",
            output_name="consensus_tree",
            artifact_kind="bootstrap-consensus-tree",
        )
        _validate_matching_tree_taxa(
            engine_name="iqtree",
            workflow="bootstrap-consensus",
            reference_tree_path=consensus_tree_path,
            comparison_tree_set_path=bootstrap_trees_path,
            reference_output_name="consensus_tree",
            comparison_output_name="bootstrap_trees",
            artifact_kind="bootstrap-tree-set",
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=None,
            support_tree_path=consensus_tree_path,
        )
        _validate_support_value_count(
            engine_name="iqtree",
            workflow="bootstrap-consensus",
            path=consensus_tree_path,
            output_name="consensus_tree",
            artifact_kind="bootstrap-consensus-tree",
            support_value_count=iqtree_summary.support_value_count,
            support_kind="bootstrap consensus support",
        )
        _validate_complete_support_coverage(
            engine_name="iqtree",
            workflow="bootstrap-consensus",
            path=consensus_tree_path,
            output_name="consensus_tree",
            artifact_kind="bootstrap-consensus-tree",
            annotated_branch_count=iqtree_summary.support_value_count,
            support_kind="bootstrap consensus support",
        )
        model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow="bootstrap-consensus",
        engine_name="iqtree",
        input_paths=[bootstrap_trees_path],
        output_paths=_existing_iqtree_outputs(prefix_path, include_consensus=True),
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([bootstrap_trees_path]),
        output_checksums={},
        config={
            "minimum_support": minimum_support,
            "timeout_seconds": timeout_seconds,
        },
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        notes=[
            "consensus tree validated as parseable Newick output",
            *(
                []
                if iqtree_summary.support_value_count == 0
                else ["support values parsed from the bootstrap consensus tree"]
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)
