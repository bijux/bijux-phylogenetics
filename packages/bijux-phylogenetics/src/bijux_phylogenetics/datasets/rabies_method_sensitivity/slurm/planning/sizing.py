from __future__ import annotations

from math import ceil
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmJobPlanRow
from .interfaces import TaskRecordLike, VariantRunLike
from .shared import (
    MEBIBYTE,
    MINIMUM_MEMORY_MIB,
    MINIMUM_SCRATCH_MIB,
    MINIMUM_WALLCLOCK_MINUTES,
    directory_bytes,
    format_slurm_time,
    round_up,
)


def build_job_plan_row(
    *,
    dataset_id: str,
    workflow_prefix: str,
    taxon_count: int,
    variant_run: VariantRunLike,
    task_record: TaskRecordLike,
    iqtree_threads: int,
    bootstrap_replicates: int,
) -> RabiesMethodSensitivitySlurmJobPlanRow:
    job_name = f"{workflow_prefix}-{variant_run.config.variant_id}"[:128]
    estimated_cpus_per_task = max(
        iqtree_threads,
        2 if variant_run.config.alignment_mode == "ginsi" else 1,
    )
    dataset_size_score = taxon_count * max(variant_run.trimmed_alignment_length, 1)
    dataset_size_class = classify_dataset_size(dataset_size_score)
    dataset_scale = max(
        1.0,
        (taxon_count * max(variant_run.trimmed_alignment_length, 1)) / 10000,
    )
    alignment_factor = 2.0 if variant_run.config.alignment_mode == "ginsi" else 1.0
    trimming_factor = 1.15 if variant_run.config.trimming_mode == "gappyout" else 1.0
    bootstrap_scale = max(1.0, bootstrap_replicates / 1000)
    estimated_wallclock_minutes = max(
        MINIMUM_WALLCLOCK_MINUTES,
        ceil(
            14
            + dataset_scale
            * bootstrap_scale
            * ((6 * alignment_factor) + (4 * trimming_factor))
        ),
    )
    memory_request = (
        768
        + ceil(dataset_scale * 160)
        + ((estimated_cpus_per_task - 1) * 256)
        + (64 if variant_run.config.trimming_mode == "gappyout" else 0)
    )
    estimated_memory_mib = max(
        MINIMUM_MEMORY_MIB, round_up(memory_request, quantum=256)
    )
    resource_class = classify_resource_class(
        estimated_cpus_per_task=estimated_cpus_per_task,
        estimated_memory_mib=estimated_memory_mib,
        estimated_wallclock_minutes=estimated_wallclock_minutes,
    )
    observed_output_bytes = directory_bytes(task_record.output_root)
    estimated_output_mib = max(1, ceil(observed_output_bytes / MEBIBYTE))
    scratch_request = (
        128 + ceil((observed_output_bytes * 64) / MEBIBYTE) + ceil(dataset_scale * 32)
    )
    estimated_scratch_mib = max(
        MINIMUM_SCRATCH_MIB, round_up(scratch_request, quantum=128)
    )
    estimated_core_hours = round(
        (estimated_cpus_per_task * estimated_wallclock_minutes) / 60, 2
    )
    task_log_path = Path(
        "parallel-logs", f"{variant_run.config.variant_id}.log"
    ).as_posix()
    return RabiesMethodSensitivitySlurmJobPlanRow(
        dataset_id=dataset_id,
        variant_id=variant_run.config.variant_id,
        job_name=job_name,
        taxon_count=taxon_count,
        dataset_size_class=dataset_size_class,
        method_group=method_group_name(variant_run.config.alignment_mode),
        resource_class=resource_class,
        alignment_mode=variant_run.config.alignment_mode,
        trimming_mode=variant_run.config.trimming_mode,
        aligned_site_count=variant_run.alignment_length,
        trimmed_site_count=variant_run.trimmed_alignment_length,
        bootstrap_replicates=bootstrap_replicates,
        iqtree_threads=iqtree_threads,
        estimated_cpus_per_task=estimated_cpus_per_task,
        estimated_memory_mib=estimated_memory_mib,
        estimated_wallclock_minutes=estimated_wallclock_minutes,
        slurm_time=format_slurm_time(estimated_wallclock_minutes),
        estimated_scratch_mib=estimated_scratch_mib,
        observed_output_bytes=observed_output_bytes,
        estimated_output_mib=estimated_output_mib,
        estimated_core_hours=estimated_core_hours,
        bundle_output_directory=Path(
            "variants", variant_run.config.variant_id
        ).as_posix(),
        task_log_path=task_log_path,
        suggested_sbatch_options=(
            f"--job-name={job_name} "
            f"--cpus-per-task={estimated_cpus_per_task} "
            f"--mem={estimated_memory_mib}M "
            f"--time={format_slurm_time(estimated_wallclock_minutes)} "
            f"--output=slurm-logs/{variant_run.config.variant_id}.%j.out "
            f"--error=slurm-logs/{variant_run.config.variant_id}.%j.err"
        ),
    )


def classify_dataset_size(size_score: int) -> str:
    if size_score <= 20_000:
        return "compact"
    if size_score <= 200_000:
        return "moderate"
    return "large"


def classify_resource_class(
    *,
    estimated_cpus_per_task: int,
    estimated_memory_mib: int,
    estimated_wallclock_minutes: int,
) -> str:
    if (
        estimated_cpus_per_task <= 1
        and estimated_memory_mib <= 1280
        and estimated_wallclock_minutes <= 30
    ):
        return "standard"
    if (
        estimated_cpus_per_task <= 2
        and estimated_memory_mib <= 2048
        and estimated_wallclock_minutes <= 45
    ):
        return "elevated"
    return "heavy"


def method_group_name(alignment_mode: str) -> str:
    return f"mafft-{alignment_mode}"
