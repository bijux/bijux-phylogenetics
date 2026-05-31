from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from bijux_phylogenetics.bayesian.posterior_execution import (
    run_bayesian_posterior_execution,
)

from .logs import (
    parse_beast_log,
)
from .posterior_trees import (
    parse_beast_posterior_tree_samples,
)

if TYPE_CHECKING:
    from bijux_phylogenetics.engines.workflows.models import EngineWorkflowReport


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
    from bijux_phylogenetics.engines.validation.preflight import (
        require_external_engine_surface,
    )

    from ._shared import (
        _beast_artifact_error,
        _beast_output_path,
        read_engine_version,
        resolve_engine_executable,
        validate_timeout_seconds,
    )

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
    require_external_engine_surface(
        workflow_id="beast-posterior",
        summary="BEAST posterior inference workflow.",
        required_engines=("beast",),
        executables={"beast": executable},
        preserve_missing_error=True,
    )
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

    def _validate_consistent_posterior_states() -> None:
        from ._shared import _beast_artifact_error

        log_report = parse_beast_log(posterior_log_path)
        tree_report = parse_beast_posterior_tree_samples(
            posterior_trees_path,
            burnin_fraction=0.0,
        )
        log_states = {row.state for row in log_report.rows}
        tree_states = set(tree_report.sampled_states)
        if tree_states.issubset(log_states):
            return
        raise _beast_artifact_error(
            f"BEAST posterior outputs disagree on sampled states: {posterior_trees_path}",
            code="beast_outputs_inconsistent_states",
            path=posterior_trees_path,
            artifact_kind="beast-posterior-trees",
            details={
                "log_path": str(posterior_log_path),
                "log_states": sorted(log_states),
                "tree_states": sorted(tree_states),
                "missing_log_states": sorted(tree_states - log_states),
            },
        )

    def validate_outputs() -> None:
        _validate_consistent_posterior_states()

    return run_bayesian_posterior_execution(
        engine_name="BEAST",
        executable=resolved,
        version=version,
        command=command,
        input_paths=[xml_path],
        output_paths={
            "posterior_log": posterior_log_path,
            "posterior_trees": posterior_trees_path,
        },
        manifest_path=manifest_path,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        work_dir=xml_path.parent,
        timeout_seconds=timeout_seconds,
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
        ],
        resume=resume,
        incomplete_run_policy=incomplete_run_policy,
        validate_outputs=validate_outputs,
    )
