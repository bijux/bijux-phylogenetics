from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import TypeAlias

from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from ..models import (
    RabiesMethodSensitivityPanelDataset,
    RabiesMethodSensitivityTaskRecord,
    RabiesMethodSensitivityVariant,
    RabiesMethodSensitivityVariantRun,
)

VariantWorkflowRunner: TypeAlias = Callable[..., RabiesMethodSensitivityVariantRun]


def _run_workflow_tasks(
    *,
    dataset: RabiesMethodSensitivityPanelDataset,
    out_dir: Path,
    execution_mode: str,
    parallel_workers: int,
    mafft_executable: str | Path,
    trimal_executable: str | Path,
    iqtree_executable: str | Path,
    fasttree_executable: str | Path,
    iqtree_seed: int,
    iqtree_threads: int,
    bootstrap_replicates: int,
    run_variant_workflow: VariantWorkflowRunner,
    write_task_log: Callable[..., None],
) -> tuple[
    tuple[RabiesMethodSensitivityTaskRecord, ...],
    dict[str, RabiesMethodSensitivityVariantRun],
]:
    task_log_root = out_dir / "parallel-logs"
    task_log_root.mkdir(parents=True, exist_ok=True)
    task_records_by_variant: dict[str, RabiesMethodSensitivityTaskRecord] = {}
    variant_runs_by_variant: dict[str, RabiesMethodSensitivityVariantRun] = {}

    def run_variant(
        variant: RabiesMethodSensitivityVariant,
    ) -> tuple[
        RabiesMethodSensitivityTaskRecord,
        RabiesMethodSensitivityVariantRun | None,
    ]:
        return _run_variant_task(
            dataset=dataset,
            variant=variant,
            out_dir=out_dir,
            execution_mode=execution_mode,
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            fasttree_executable=fasttree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
            run_variant_workflow=run_variant_workflow,
            write_task_log=write_task_log,
            task_log_root=task_log_root,
        )

    if execution_mode == "serial":
        for variant in dataset.variants:
            task_record, variant_run = run_variant(variant)
            task_records_by_variant[variant.variant_id] = task_record
            if variant_run is not None:
                variant_runs_by_variant[variant.variant_id] = variant_run
    else:
        with ThreadPoolExecutor(max_workers=parallel_workers) as executor:
            futures = {
                executor.submit(run_variant, variant): variant.variant_id
                for variant in dataset.variants
            }
            for future in as_completed(futures):
                task_record, variant_run = future.result()
                task_records_by_variant[task_record.variant_id] = task_record
                if variant_run is not None:
                    variant_runs_by_variant[task_record.variant_id] = variant_run

    task_records = tuple(
        task_records_by_variant[variant.variant_id] for variant in dataset.variants
    )
    return task_records, variant_runs_by_variant


def _run_variant_task(
    *,
    dataset: RabiesMethodSensitivityPanelDataset,
    variant: RabiesMethodSensitivityVariant,
    out_dir: Path,
    execution_mode: str,
    mafft_executable: str | Path,
    trimal_executable: str | Path,
    iqtree_executable: str | Path,
    fasttree_executable: str | Path,
    iqtree_seed: int,
    iqtree_threads: int,
    bootstrap_replicates: int,
    run_variant_workflow: VariantWorkflowRunner,
    write_task_log: Callable[..., None],
    task_log_root: Path,
) -> tuple[
    RabiesMethodSensitivityTaskRecord,
    RabiesMethodSensitivityVariantRun | None,
]:
    variant_root = out_dir / "variants" / variant.variant_id
    variant_log_path = task_log_root / f"{variant.variant_id}.log"
    try:
        variant_run = run_variant_workflow(
            dataset=dataset,
            variant=variant,
            variant_root=variant_root,
            mafft_executable=mafft_executable,
            trimal_executable=trimal_executable,
            iqtree_executable=iqtree_executable,
            fasttree_executable=fasttree_executable,
            iqtree_seed=iqtree_seed,
            iqtree_threads=iqtree_threads,
            bootstrap_replicates=bootstrap_replicates,
        )
    except Exception as error:
        error_code = (
            error.code
            if isinstance(error, PhylogeneticsError) and error.code is not None
            else "parallel_variant_failed"
        )
        error_message = str(error)
        write_task_log(
            variant_log_path,
            variant=variant,
            execution_mode=execution_mode,
            status="failed",
            output_root=Path("variants") / variant.variant_id,
            error_code=error_code,
            error_message=error_message,
        )
        return (
            RabiesMethodSensitivityTaskRecord(
                variant_id=variant.variant_id,
                label=variant.label,
                execution_mode=execution_mode,
                status="failed",
                log_path=variant_log_path,
                output_root=variant_root,
                error_code=error_code,
                error_message=error_message,
            ),
            None,
        )
    write_task_log(
        variant_log_path,
        variant=variant,
        execution_mode=execution_mode,
        status="succeeded",
        output_root=Path("variants") / variant.variant_id,
        error_code=None,
        error_message=None,
    )
    return (
        RabiesMethodSensitivityTaskRecord(
            variant_id=variant.variant_id,
            label=variant.label,
            execution_mode=execution_mode,
            status="succeeded",
            log_path=variant_log_path,
            output_root=variant_root,
            error_code=None,
            error_message=None,
        ),
        variant_run,
    )
