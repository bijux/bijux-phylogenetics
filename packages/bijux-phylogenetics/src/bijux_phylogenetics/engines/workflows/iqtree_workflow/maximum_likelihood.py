from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet
from bijux_phylogenetics.runtime.errors import EngineWorkflowError, PhylogeneticsError

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
    _ensure_inference_ready_alignment,
    _persist_workflow_report,
    _prefix_path,
    _record_output_validation_failure,
    _resolve_incomplete_workflow_state,
    _resume_existing_workflow,
    _validate_iqtree_required_artifacts,
    _validate_tree_output,
)
from .partitions import _prepare_iqtree_partitions
from .shared import (
    _build_iqtree_model_selection_summary,
    _build_iqtree_summary,
    _existing_iqtree_outputs,
    _iqtree_execution_controls,
    _iqtree_partition_supports_fixed_model,
    _iqtree_sequence_type_flag,
    _validate_iqtree_model_result,
)


def run_maximum_likelihood_tree_inference(
    input_path: Path,
    *,
    out_dir: Path,
    model: str,
    prefix: str = "maximum-likelihood",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run an external maximum-likelihood tree inference workflow."""
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    require_external_engine_surface(
        workflow_id="maximum-likelihood-tree",
        summary="IQ-TREE maximum-likelihood tree-inference workflow.",
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
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
    if prepared_partitions is not None and not _iqtree_partition_supports_fixed_model(
        model=model,
        mixed_data_types=prepared_partitions.mixed_data_types,
    ):
        raise EngineWorkflowError(
            "mixed DNA/protein partition analyses require a model-selection keyword such as MF, MFP, TEST, or TESTMERGE"
        )
    tree_path = prefix_path.with_suffix(".treefile")
    report_path = prefix_path.with_suffix(".iqtree")
    log_path = prefix_path.with_suffix(".log")
    command = [
        resolved,
        *(
            prepared_partitions.command_args
            if prepared_partitions is not None
            else ["-s", str(input_path.resolve())]
        ),
        *(
            []
            if prepared_partitions is not None and prepared_partitions.mixed_data_types
            else _iqtree_sequence_type_flag(input_path, sequence_type)
        ),
        *_iqtree_execution_controls(seed=seed, threads=threads),
        "-m",
        model,
        "-pre",
        str(prefix_path.resolve()),
    ]
    if resume:
        resumed = _resume_existing_workflow(
            manifest_path=manifest_path,
            input_paths=(
                [input_path] if partition_path is None else [input_path, partition_path]
            ),
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
        workflow="maximum-likelihood-tree",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "tree": tree_path,
            "iqtree_report": report_path,
            "iqtree_log": log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_iqtree_required_artifacts(
            prefix_path,
            workflow="maximum-likelihood-tree",
        )
        _validate_tree_output(
            tree_path,
            engine_name="iqtree",
            workflow="maximum-likelihood-tree",
            output_name="tree",
            artifact_kind="maximum-likelihood-tree",
        )
        selected_model = _validate_iqtree_model_result(
            prefix_path,
            workflow="maximum-likelihood-tree",
            default_selected_model=model,
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=selected_model,
            support_tree_path=tree_path,
        )
        model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    report = EngineWorkflowReport(
        workflow="maximum-likelihood-tree",
        engine_name="iqtree",
        input_paths=(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_paths={
            **({} if prepared_partitions is None else prepared_partitions.output_paths),
            **_existing_iqtree_outputs(prefix_path, include_tree=True),
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        config={
            "model": model,
            "sequence_type": sequence_type,
            "partition_path": None if partition_path is None else str(partition_path),
            "seed": seed,
            "threads": threads,
            "timeout_seconds": timeout_seconds,
        },
        selected_model=iqtree_summary.selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        notes=[
            *([] if prepared_partitions is None else prepared_partitions.notes),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            "maximum-likelihood tree validated as parseable Newick output",
            *(
                []
                if iqtree_summary.log_likelihood is None
                else ["log-likelihood parsed from iqtree inference artifacts"]
            ),
            *(
                []
                if iqtree_summary.support_value_count == 0
                else ["support values parsed from the inferred maximum-likelihood tree"]
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)
