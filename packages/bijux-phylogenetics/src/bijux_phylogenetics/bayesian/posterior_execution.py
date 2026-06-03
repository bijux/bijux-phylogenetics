from __future__ import annotations

from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any

from bijux_phylogenetics.runtime.errors import PhylogeneticsError

if TYPE_CHECKING:
    from bijux_phylogenetics.engines.common import EngineVersionInfo
    from bijux_phylogenetics.engines.workflows.models import EngineWorkflowReport

PosteriorOutputValidator = Callable[[], None]


def run_bayesian_posterior_execution(
    *,
    engine_name: str,
    executable: str,
    version: EngineVersionInfo,
    command: list[str],
    input_paths: list[Path],
    output_paths: dict[str, Path],
    manifest_path: Path,
    stdout_path: Path,
    stderr_path: Path,
    work_dir: Path,
    timeout_seconds: float | None,
    config: dict[str, Any],
    notes: list[str],
    resume: bool,
    incomplete_run_policy: str,
    validate_outputs: PosteriorOutputValidator,
    workflow: str = "posterior-tree-inference",
) -> EngineWorkflowReport:
    """Run one Bayesian posterior workflow with identical resume safety rules."""
    from bijux_phylogenetics.engines.common import (
        build_file_checksums,
        execute_engine_command,
    )
    from bijux_phylogenetics.engines.workflows.models import EngineWorkflowReport
    from bijux_phylogenetics.engines.workflows.state import (
        _persist_workflow_report,
        _record_output_validation_failure,
        _resolve_incomplete_workflow_state,
        _resume_existing_workflow,
    )

    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=input_paths,
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
        engine_name=engine_name,
        workflow=workflow,
        executable=executable,
        version=version,
        command_args=command[1:],
        work_dir=work_dir,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        output_paths=output_paths,
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        validate_outputs()
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow=workflow,
        engine_name=engine_name,
        input_paths=input_paths,
        output_paths=output_paths,
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(input_paths),
        output_checksums={},
        config=config,
        notes=[*notes, *incomplete_notes],
    )
    return _persist_workflow_report(report)
