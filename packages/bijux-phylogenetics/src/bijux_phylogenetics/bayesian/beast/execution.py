from __future__ import annotations

from ._shared import (
    EngineWorkflowReport,
    Path,
    PhylogeneticsError,
    _beast_artifact_error,
    _beast_output_path,
    _persist_workflow_report,
    _record_output_validation_failure,
    _resolve_incomplete_workflow_state,
    _resume_existing_workflow,
    build_file_checksums,
    execute_engine_command,
    read_engine_version,
    resolve_engine_executable,
    validate_timeout_seconds,
)
from .logs import (
    parse_beast_log,
)
from .posterior_trees import (
    parse_beast_posterior_tree_samples,
)


def run_beast_posterior_inference(
    xml_path: Path,
    *,
    executable: str | Path = "beast",
    overwrite: bool = True,
    threads: int = 1,
    seed: int = 1,
    resume: bool = False,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run a prepared BEAST XML analysis and validate the primary posterior outputs."""
    if not xml_path.exists():
        raise _beast_artifact_error(
            f"BEAST analysis XML file was not found: {xml_path}",
            code="beast_xml_missing_file",
            path=xml_path,
            artifact_kind="beast-analysis-xml",
        )
    validate_timeout_seconds(timeout_seconds)
    if threads < 1:
        raise ValueError(f"threads must be positive, got {threads}")
    if seed < 1:
        raise ValueError(f"seed must be positive, got {seed}")
    resolved = resolve_engine_executable(executable)
    manifest_path = xml_path.with_suffix(".manifest.json")
    stdout_path = xml_path.with_suffix(".stdout.log")
    stderr_path = xml_path.with_suffix(".stderr.log")
    posterior_log_path = _beast_output_path(xml_path, seed=seed, suffix="log")
    posterior_trees_path = _beast_output_path(xml_path, seed=seed, suffix="trees")
    version = read_engine_version(
        "BEAST",
        executable,
        version_args=("-version",),
        timeout_seconds=timeout_seconds,
    )
    command = [
        resolved,
        *(["-overwrite"] if overwrite else []),
        "-threads",
        str(threads),
        "-seed",
        str(seed),
        xml_path.name,
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=[xml_path],
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
        engine_name="BEAST",
        workflow="posterior-tree-inference",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=xml_path.parent,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        output_paths={
            "posterior_log": posterior_log_path,
            "posterior_trees": posterior_trees_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        parse_beast_log(posterior_log_path)
        parse_beast_posterior_tree_samples(posterior_trees_path, burnin_fraction=0.0)
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow="posterior-tree-inference",
        engine_name="BEAST",
        input_paths=[xml_path],
        output_paths={
            "posterior_log": posterior_log_path,
            "posterior_trees": posterior_trees_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums([xml_path]),
        output_checksums={},
        config={
            "threads": threads,
            "seed": seed,
            "overwrite": overwrite,
            "timeout_seconds": timeout_seconds,
        },
        notes=[
            "BEAST posterior log and posterior tree set validated after engine execution",
            f"beast threads: {threads}",
            f"beast random seed: {seed}",
            *(
                ["existing posterior outputs are overwritten before engine execution"]
                if overwrite
                else []
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)
