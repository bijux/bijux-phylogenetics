from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from ...artifacts.iqtree import write_iqtree_model_candidates_table
from ...common import (
    build_engine_output_error,
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
)
from .partitions import _prepare_iqtree_partitions
from .shared import (
    _build_iqtree_model_selection_summary,
    _build_iqtree_summary,
    _existing_iqtree_outputs,
    _iqtree_execution_controls,
    _iqtree_sequence_type_flag,
    _validate_iqtree_model_result,
)


def run_model_selection(
    input_path: Path,
    *,
    out_dir: Path,
    prefix: str = "model-selection",
    executable: str | Path = "iqtree2",
    sequence_type: AlignmentAlphabet | None = None,
    partition_path: Path | None = None,
    resume: bool = False,
    seed: int = 1,
    threads: int = 1,
    timeout_seconds: float | None = None,
    incomplete_run_policy: str = "reject",
) -> EngineWorkflowReport:
    """Run a model-selection workflow on an aligned FASTA file."""
    _ensure_inference_ready_alignment(input_path)
    validate_timeout_seconds(timeout_seconds)
    require_external_engine_surface(
        workflow_id="model-selection",
        summary="IQ-TREE model-selection workflow.",
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
    iqtree_report_path = prefix_path.with_suffix(".iqtree")
    iqtree_log_path = prefix_path.with_suffix(".log")
    prepared_partitions = (
        None
        if partition_path is None
        else _prepare_iqtree_partitions(
            input_path,
            partition_path,
            prefix_path=prefix_path,
        )
    )
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
        "MF",
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
        workflow="model-selection",
        executable=resolved,
        version=version,
        command_args=command[1:],
        work_dir=out_dir,
        stdout_path=prefix_path.with_suffix(".stdout.log"),
        stderr_path=prefix_path.with_suffix(".stderr.log"),
        output_paths={
            "iqtree_report": iqtree_report_path,
            "iqtree_log": iqtree_log_path,
        },
        manifest_path=manifest_path,
        timeout_seconds=timeout_seconds,
    )
    try:
        _validate_iqtree_required_artifacts(prefix_path, workflow="model-selection")
        selected_model = _validate_iqtree_model_result(
            prefix_path,
            workflow="model-selection",
        )
        iqtree_summary = _build_iqtree_summary(
            prefix_path,
            default_selected_model=selected_model,
        )
        model_selection_summary = _build_iqtree_model_selection_summary(prefix_path)
        if (
            model_selection_summary is None
            or model_selection_summary.candidate_count < 1
        ):
            raise build_engine_output_error(
                "iqtree model-selection did not expose a parsable candidate-model table",
                code="iqtree_model_candidates_missing",
                engine_name="iqtree",
                workflow="model-selection",
                path=iqtree_report_path,
                output_name="iqtree_report",
                artifact_kind="iqtree-model-candidates",
            )
    except PhylogeneticsError as error:
        _record_output_validation_failure(manifest_path, run, error)
        raise
    selected_model_path = prefix_path.with_suffix(".selected-model.txt")
    selected_model_path.write_text(selected_model + "\n", encoding="utf-8")
    model_candidates_path = prefix_path.with_suffix(".model-candidates.tsv")
    write_iqtree_model_candidates_table(model_candidates_path, model_selection_summary)
    report = EngineWorkflowReport(
        workflow="model-selection",
        engine_name="iqtree",
        input_paths=(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_paths={
            **({} if prepared_partitions is None else prepared_partitions.output_paths),
            **_existing_iqtree_outputs(prefix_path),
            "selected_model": selected_model_path,
            "model_candidates": model_candidates_path,
        },
        run=run,
        manifest_path=manifest_path,
        input_checksums=build_file_checksums(
            [input_path] if partition_path is None else [input_path, partition_path]
        ),
        output_checksums={},
        config={
            "sequence_type": sequence_type,
            "partition_path": None if partition_path is None else str(partition_path),
            "seed": seed,
            "threads": threads,
            "timeout_seconds": timeout_seconds,
        },
        selected_model=selected_model,
        log_likelihood=iqtree_summary.log_likelihood,
        iqtree_summary=iqtree_summary,
        model_selection_summary=model_selection_summary,
        notes=[
            *([] if prepared_partitions is None else prepared_partitions.notes),
            f"iqtree random seed: {seed}",
            f"iqtree threads: {threads}",
            "best-fit substitution model parsed from engine output",
            f"parsed {model_selection_summary.candidate_count} candidate substitution models from iqtree output",
            *(
                []
                if model_selection_summary.selected_criterion is None
                else [
                    "model-selection workflow exposed the governing information criterion"
                ]
            ),
            *(
                []
                if iqtree_summary.log_likelihood is None
                else [
                    "model-selection workflow exposed a parsable log-likelihood score"
                ]
            ),
            *incomplete_notes,
        ],
    )
    return _persist_workflow_report(report)
