from __future__ import annotations

import os
from pathlib import Path

from bijux_phylogenetics.engines.common import (
    EngineActiveRunRecord,
    acquire_active_engine_run,
    clear_active_engine_run,
    utc_now_text,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

from ..config import (
    _resolve_selected_variant_dataset,
    load_rabies_method_sensitivity_panel_dataset,
)
from ..models import (
    RabiesMethodSensitivityPanelWorkflowReport,
)
from .clade_review import _aggregate_clades
from .conclusions import _build_conclusion_rows
from .records import (
    _format_float,
    _workflow_execution_record_path,
    _write_task_log,
    _write_workflow_execution_record,
)
from .task_execution import _run_workflow_tasks
from .topology_review import _build_preprocessing_comparison_rows
from .variant_execution import _run_variant_workflow

__all__ = ["run_rabies_method_sensitivity_panel_workflow"]


def run_rabies_method_sensitivity_panel_workflow(
    out_dir: Path,
    *,
    mafft_executable: str | Path = "mafft",
    trimal_executable: str | Path = "trimal",
    iqtree_executable: str | Path = "iqtree2",
    fasttree_executable: str | Path = "FastTree",
    iqtree_seed: int | None = None,
    iqtree_threads: int | None = None,
    bootstrap_replicates: int | None = None,
    parallel_workers: int | None = None,
    variant_ids: tuple[str, ...] | None = None,
) -> RabiesMethodSensitivityPanelWorkflowReport:
    """Run the owned method-sensitivity workflow over the packaged rabies panel."""
    dataset = _resolve_selected_variant_dataset(
        load_rabies_method_sensitivity_panel_dataset(),
        variant_ids=variant_ids,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    started_at_utc = utc_now_text()
    resolved_seed = dataset.iqtree_seed if iqtree_seed is None else iqtree_seed
    resolved_threads = (
        dataset.iqtree_threads if iqtree_threads is None else iqtree_threads
    )
    resolved_bootstrap_replicates = (
        dataset.bootstrap_replicates
        if bootstrap_replicates is None
        else bootstrap_replicates
    )
    resolved_parallel_workers = (
        dataset.parallel_workers if parallel_workers is None else parallel_workers
    )
    if resolved_parallel_workers < 1:
        raise ValueError(
            f"parallel_workers must be at least 1, got {resolved_parallel_workers}"
        )
    execution_mode = (
        "parallel"
        if resolved_parallel_workers > 1 and len(dataset.variants) > 1
        else "serial"
    )
    execution_record_path = _workflow_execution_record_path(out_dir)
    active_record = EngineActiveRunRecord(
        engine_name="bijux",
        workflow="rabies_method_sensitivity_panel",
        executable="bijux-phylogenetics",
        working_directory=out_dir,
        manifest_path=execution_record_path,
        command=["run_rabies_method_sensitivity_panel_workflow"],
        process_id=os.getpid(),
        started_at_utc=started_at_utc,
    )
    try:
        acquire_active_engine_run(active_record)
    except EngineWorkflowError as error:
        if error.code == "engine_workflow_already_running":
            raise EngineWorkflowError(
                "rabies method-sensitivity workflow is already running for the requested output root",
                code="rabies_method_sensitivity_workflow_already_running",
                details={
                    "output_root": str(out_dir),
                    "execution_record_path": str(execution_record_path),
                    **(error.details or {}),
                },
            ) from error
        if error.code == "engine_workflow_running_marker_invalid":
            raise EngineWorkflowError(
                "rabies method-sensitivity workflow could not verify the active-run marker for the requested output root",
                code="rabies_method_sensitivity_workflow_running_marker_invalid",
                details={
                    "output_root": str(out_dir),
                    "execution_record_path": str(execution_record_path),
                    **(error.details or {}),
                },
            ) from error
        raise

    try:
        task_records, variant_runs_by_variant = _run_workflow_tasks(
            dataset=dataset,
            out_dir=out_dir,
            execution_mode=execution_mode,
            parallel_workers=resolved_parallel_workers,
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            fasttree_executable=fasttree_executable,
            iqtree_seed=resolved_seed,
            iqtree_threads=resolved_threads,
            bootstrap_replicates=resolved_bootstrap_replicates,
            run_variant_workflow=_run_variant_workflow,
            write_task_log=_write_task_log,
        )
        failed_task_records = [
            record for record in task_records if record.status != "succeeded"
        ]
        ended_at_utc = utc_now_text()
        if failed_task_records:
            _write_workflow_execution_record(
                execution_record_path,
                dataset=dataset,
                execution_mode=execution_mode,
                parallel_workers=resolved_parallel_workers,
                task_records=task_records,
                started_at_utc=started_at_utc,
                ended_at_utc=ended_at_utc,
                status="failed",
            )
            raise EngineWorkflowError(
                "rabies method-sensitivity workflow left one or more parallel tasks failed while preserving successful isolated outputs",
                code="workflow_parallel_task_failed",
                details={
                    "failed_variants": [
                        record.variant_id for record in failed_task_records
                    ],
                    "successful_variants": [
                        record.variant_id
                        for record in task_records
                        if record.status == "succeeded"
                    ],
                    "parallel_workers": resolved_parallel_workers,
                    "execution_mode": execution_mode,
                    "task_logs": {
                        record.variant_id: str(record.log_path)
                        for record in task_records
                    },
                    "execution_record_path": str(execution_record_path),
                },
            )
        variant_runs = [
            variant_runs_by_variant[variant.variant_id] for variant in dataset.variants
        ]

        preprocessing_comparison_rows = tuple(
            _build_preprocessing_comparison_rows(variant_runs)
        )
        stable_clade_rows = tuple(_aggregate_clades(variant_runs, stable_only=True))
        changed_clade_rows = tuple(_aggregate_clades(variant_runs, stable_only=False))
        conclusion_rows = tuple(
            _build_conclusion_rows(
                dataset=dataset,
                variant_runs=variant_runs,
                preprocessing_comparison_rows=preprocessing_comparison_rows,
            )
        )
        report = RabiesMethodSensitivityPanelWorkflowReport(
            dataset=dataset,
            execution_record_path=execution_record_path,
            iqtree_seed=resolved_seed,
            iqtree_threads=resolved_threads,
            bootstrap_replicates=resolved_bootstrap_replicates,
            parallel_workers=resolved_parallel_workers,
            execution_mode=execution_mode,
            task_records=task_records,
            variant_runs=tuple(variant_runs),
            preprocessing_comparison_rows=preprocessing_comparison_rows,
            stable_clade_rows=stable_clade_rows,
            changed_clade_rows=changed_clade_rows,
            conclusion_rows=conclusion_rows,
        )
        _write_workflow_execution_record(
            execution_record_path,
            dataset=dataset,
            execution_mode=execution_mode,
            parallel_workers=resolved_parallel_workers,
            task_records=task_records,
            started_at_utc=started_at_utc,
            ended_at_utc=ended_at_utc,
            status="succeeded",
        )
        return report
    finally:
        clear_active_engine_run(execution_record_path)
