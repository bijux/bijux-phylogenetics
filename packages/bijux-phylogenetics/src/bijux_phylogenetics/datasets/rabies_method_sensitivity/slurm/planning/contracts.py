from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmAssumptionRow:
    """One explicit planning assumption used to size a Slurm job."""

    assumption_id: str
    parameter: str
    value: str
    rationale: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmJobPlanRow:
    """One Slurm-ready estimate for one governed rabies workflow variant."""

    dataset_id: str
    variant_id: str
    job_name: str
    taxon_count: int
    dataset_size_class: str
    method_group: str
    resource_class: str
    alignment_mode: str
    trimming_mode: str
    aligned_site_count: int
    trimmed_site_count: int
    bootstrap_replicates: int
    iqtree_threads: int
    estimated_cpus_per_task: int
    estimated_memory_mib: int
    estimated_wallclock_minutes: int
    slurm_time: str
    estimated_scratch_mib: int
    observed_output_bytes: int
    estimated_output_mib: int
    estimated_core_hours: float
    bundle_output_directory: str
    task_log_path: str
    suggested_sbatch_options: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmPlanningReport:
    """One reviewer-facing Slurm job-planning report for the rabies batch workflow."""

    dataset_id: str
    workflow_prefix: str
    job_count: int
    bootstrap_replicates: int
    iqtree_threads: int
    total_estimated_core_hours: float
    maximum_estimated_memory_mib: int
    maximum_estimated_wallclock_minutes: int
    total_estimated_scratch_mib: int
    total_estimated_output_mib: int
    assumptions: tuple[RabiesMethodSensitivitySlurmAssumptionRow, ...]
    rows: tuple[RabiesMethodSensitivitySlurmJobPlanRow, ...]
