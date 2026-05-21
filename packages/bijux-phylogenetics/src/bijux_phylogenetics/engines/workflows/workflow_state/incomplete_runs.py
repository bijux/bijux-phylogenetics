from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError

from ...common import (
    EngineRunReport,
    active_engine_run_is_live,
    cleanup_incomplete_engine_run,
    engine_active_marker_path,
    engine_incomplete_marker_path,
    load_active_engine_run,
    load_incomplete_engine_run,
    observe_engine_outputs,
    update_incomplete_engine_run,
)

_INCOMPLETE_RUN_POLICIES = {"reject", "clean"}


def _validate_incomplete_run_policy(policy: str) -> str:
    if policy not in _INCOMPLETE_RUN_POLICIES:
        available = ", ".join(sorted(_INCOMPLETE_RUN_POLICIES))
        raise ValueError(
            f"incomplete_run_policy must be one of: {available}; got {policy}"
        )
    return policy


def _resolve_incomplete_workflow_state(
    *,
    manifest_path: Path,
    incomplete_run_policy: str,
) -> list[str]:
    _validate_incomplete_run_policy(incomplete_run_policy)
    active_record = load_active_engine_run(manifest_path)
    if active_record is not None and active_engine_run_is_live(active_record):
        raise EngineWorkflowError(
            "engine workflow is already running for the requested output manifest",
            code="engine_workflow_already_running",
            details={
                "manifest_path": str(manifest_path),
                "marker_path": str(engine_active_marker_path(manifest_path)),
                "running_process_id": active_record.process_id,
                "running_workflow": active_record.workflow,
                "running_engine_name": active_record.engine_name,
            },
        )
    record = load_incomplete_engine_run(manifest_path)
    if record is None:
        return []
    if incomplete_run_policy == "clean":
        cleanup_incomplete_engine_run(manifest_path)
        return [
            "removed outputs from a previously incomplete engine run before restarting"
        ]
    marker_path = engine_incomplete_marker_path(manifest_path)
    observed_outputs = [
        {
            "output_name": observation.output_name,
            "path": str(observation.path),
            "exists": observation.exists,
            "path_kind": observation.path_kind,
            "size_bytes": observation.size_bytes,
            "sha256": observation.sha256,
        }
        for observation in record.observed_outputs
    ]
    raise EngineWorkflowError(
        "a previous engine run left incomplete outputs and resume could not safely "
        f"reuse them; marker: {marker_path}",
        code="engine_incomplete_outputs_present",
        details={
            "manifest_path": str(manifest_path),
            "marker_path": str(marker_path),
            "engine_name": record.engine_name,
            "workflow": record.workflow,
            "failure_reason": (
                record.failure_reason
                if record.failure_reason is not None
                else "engine_run_incomplete"
            ),
            "failure_message": record.failure_message,
            "timed_out": record.timed_out,
            "exit_code": record.exit_code,
            "timeout_seconds": record.timeout_seconds,
            "missing_output_names": list(record.missing_output_names),
            "observed_outputs": observed_outputs,
            "incomplete_run_policy": incomplete_run_policy,
            "available_actions": ["resume", "clean"],
        },
    )


def _record_output_validation_failure(
    manifest_path: Path,
    run: EngineRunReport,
    error: PhylogeneticsError,
) -> None:
    observations = observe_engine_outputs(run.output_paths)
    update_incomplete_engine_run(
        manifest_path,
        ended_at_utc=run.ended_at_utc,
        timed_out=run.timed_out,
        exit_code=run.exit_code,
        failure_reason=error.code,
        failure_message=(
            f"{run.engine_name} {run.workflow} produced outputs that failed "
            f"validation: {error.code}"
        ),
        missing_output_names=[
            observation.output_name
            for observation in observations
            if not observation.exists
        ],
        observed_outputs=observations,
    )
