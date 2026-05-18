from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from math import ceil
from pathlib import Path
from typing import Protocol

__all__ = [
    "RabiesMethodSensitivitySlurmAssumptionRow",
    "RabiesMethodSensitivitySlurmJobPlanRow",
    "RabiesMethodSensitivitySlurmPlanningReport",
    "build_rabies_method_sensitivity_slurm_planning_report",
    "write_rabies_method_sensitivity_slurm_assumptions_table",
    "write_rabies_method_sensitivity_slurm_job_plan_table",
    "write_rabies_method_sensitivity_slurm_summary_json",
]

_MEBIBYTE = 1024 * 1024
_MINIMUM_MEMORY_MIB = 1024
_MINIMUM_SCRATCH_MIB = 256
_MINIMUM_WALLCLOCK_MINUTES = 20


class _VariantConfigLike(Protocol):
    variant_id: str
    alignment_mode: str
    trimming_mode: str


class _VariantRunLike(Protocol):
    config: _VariantConfigLike
    alignment_length: int
    trimmed_alignment_length: int


class _TaskRecordLike(Protocol):
    variant_id: str
    output_root: Path


class _DatasetLike(Protocol):
    dataset_id: str
    workflow_prefix: str
    taxon_count: int


class _WorkflowReportLike(Protocol):
    dataset: _DatasetLike
    task_records: tuple[_TaskRecordLike, ...]
    variant_runs: tuple[_VariantRunLike, ...]
    iqtree_threads: int
    bootstrap_replicates: int


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


def build_rabies_method_sensitivity_slurm_planning_report(
    report: _WorkflowReportLike,
) -> RabiesMethodSensitivitySlurmPlanningReport:
    """Estimate one schedulable Slurm job per declared workflow variant."""
    task_records = {record.variant_id: record for record in report.task_records}
    rows = tuple(
        _build_job_plan_row(
            dataset_id=report.dataset.dataset_id,
            workflow_prefix=report.dataset.workflow_prefix,
            taxon_count=report.dataset.taxon_count,
            variant_run=variant_run,
            task_record=task_records[variant_run.config.variant_id],
            iqtree_threads=report.iqtree_threads,
            bootstrap_replicates=report.bootstrap_replicates,
        )
        for variant_run in report.variant_runs
    )
    assumptions = (
        RabiesMethodSensitivitySlurmAssumptionRow(
            assumption_id="observed-output-footprint",
            parameter="observed_output_bytes",
            value="sum of files currently written under each variant output root",
            rationale=(
                "Output-size and scratch estimates are anchored to the real bytes "
                "written by the current governed workflow run instead of a synthetic "
                "placeholder volume."
            ),
        ),
        RabiesMethodSensitivitySlurmAssumptionRow(
            assumption_id="alignment-mode-cpu-bump",
            parameter="estimated_cpus_per_task",
            value="ginsi variants reserve at least 2 CPUs; auto variants follow iqtree_threads",
            rationale=(
                "The ginsi alignment mode is materially heavier than the auto mode, "
                "so the planner allocates an extra CPU even when the governed IQ-TREE "
                "thread count stays at 1."
            ),
        ),
        RabiesMethodSensitivitySlurmAssumptionRow(
            assumption_id="bootstrap-driven-wallclock",
            parameter="estimated_wallclock_minutes",
            value="20-minute floor plus linear scaling with trimmed sites, taxa, and bootstrap count",
            rationale=(
                "The batch wallclock estimate should stay conservative when the "
                "bootstrap burden rises, even on this compact rabies panel."
            ),
        ),
        RabiesMethodSensitivitySlurmAssumptionRow(
            assumption_id="minimum-memory-floor",
            parameter="estimated_memory_mib",
            value="1024 MiB minimum, rounded to 256 MiB blocks",
            rationale=(
                "The dataset is small, but Slurm requests still need enough room for "
                "alignment, inference, logging, and report materialization overhead."
            ),
        ),
        RabiesMethodSensitivitySlurmAssumptionRow(
            assumption_id="scratch-buffer",
            parameter="estimated_scratch_mib",
            value="observed output footprint inflated 64x with a 256 MiB minimum",
            rationale=(
                "Scratch space should cover intermediate alignments, trees, manifests, "
                "and logs without pretending that the final output footprint is the "
                "whole temporary working set."
            ),
        ),
    )
    return RabiesMethodSensitivitySlurmPlanningReport(
        dataset_id=report.dataset.dataset_id,
        workflow_prefix=report.dataset.workflow_prefix,
        job_count=len(rows),
        bootstrap_replicates=report.bootstrap_replicates,
        iqtree_threads=report.iqtree_threads,
        total_estimated_core_hours=round(
            sum(row.estimated_core_hours for row in rows), 2
        ),
        maximum_estimated_memory_mib=max(
            row.estimated_memory_mib for row in rows
        ),
        maximum_estimated_wallclock_minutes=max(
            row.estimated_wallclock_minutes for row in rows
        ),
        total_estimated_scratch_mib=sum(row.estimated_scratch_mib for row in rows),
        total_estimated_output_mib=sum(row.estimated_output_mib for row in rows),
        assumptions=assumptions,
        rows=rows,
    )


def write_rabies_method_sensitivity_slurm_job_plan_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmPlanningReport,
) -> Path:
    """Write the per-job Slurm plan as a reviewer-facing TSV table."""
    rows = [
        [
            "dataset_id",
            "variant_id",
            "job_name",
            "alignment_mode",
            "trimming_mode",
            "aligned_site_count",
            "trimmed_site_count",
            "bootstrap_replicates",
            "iqtree_threads",
            "estimated_cpus_per_task",
            "estimated_memory_mib",
            "estimated_wallclock_minutes",
            "slurm_time",
            "estimated_scratch_mib",
            "observed_output_bytes",
            "estimated_output_mib",
            "estimated_core_hours",
            "bundle_output_directory",
            "task_log_path",
            "suggested_sbatch_options",
        ]
    ]
    for row in report.rows:
        rows.append(
            [
                row.dataset_id,
                row.variant_id,
                row.job_name,
                row.alignment_mode,
                row.trimming_mode,
                str(row.aligned_site_count),
                str(row.trimmed_site_count),
                str(row.bootstrap_replicates),
                str(row.iqtree_threads),
                str(row.estimated_cpus_per_task),
                str(row.estimated_memory_mib),
                str(row.estimated_wallclock_minutes),
                row.slurm_time,
                str(row.estimated_scratch_mib),
                str(row.observed_output_bytes),
                str(row.estimated_output_mib),
                _format_float(row.estimated_core_hours),
                row.bundle_output_directory,
                row.task_log_path,
                row.suggested_sbatch_options,
            ]
        )
    return _write_tsv(path, rows)


def write_rabies_method_sensitivity_slurm_assumptions_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmPlanningReport,
) -> Path:
    """Write the sizing assumptions behind the Slurm planning report."""
    rows = [["assumption_id", "parameter", "value", "rationale"]]
    for row in report.assumptions:
        rows.append([row.assumption_id, row.parameter, row.value, row.rationale])
    return _write_tsv(path, rows)


def write_rabies_method_sensitivity_slurm_summary_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmPlanningReport,
) -> Path:
    """Write the structured Slurm planning summary for machine review."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "dataset_id": report.dataset_id,
        "workflow_prefix": report.workflow_prefix,
        "job_count": report.job_count,
        "bootstrap_replicates": report.bootstrap_replicates,
        "iqtree_threads": report.iqtree_threads,
        "total_estimated_core_hours": report.total_estimated_core_hours,
        "maximum_estimated_memory_mib": report.maximum_estimated_memory_mib,
        "maximum_estimated_wallclock_minutes": report.maximum_estimated_wallclock_minutes,
        "total_estimated_scratch_mib": report.total_estimated_scratch_mib,
        "total_estimated_output_mib": report.total_estimated_output_mib,
        "assumptions": [asdict(row) for row in report.assumptions],
        "jobs": [asdict(row) for row in report.rows],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _build_job_plan_row(
    *,
    dataset_id: str,
    workflow_prefix: str,
    taxon_count: int,
    variant_run: _VariantRunLike,
    task_record: _TaskRecordLike,
    iqtree_threads: int,
    bootstrap_replicates: int,
) -> RabiesMethodSensitivitySlurmJobPlanRow:
    job_name = (
        f"{workflow_prefix}-{variant_run.config.variant_id}"
    )[:128]
    estimated_cpus_per_task = max(
        iqtree_threads,
        2 if variant_run.config.alignment_mode == "ginsi" else 1,
    )
    dataset_scale = max(
        1.0,
        (taxon_count * max(variant_run.trimmed_alignment_length, 1)) / 10000,
    )
    alignment_factor = 2.0 if variant_run.config.alignment_mode == "ginsi" else 1.0
    trimming_factor = (
        1.15 if variant_run.config.trimming_mode == "gappyout" else 1.0
    )
    bootstrap_scale = max(1.0, bootstrap_replicates / 1000)
    estimated_wallclock_minutes = max(
        _MINIMUM_WALLCLOCK_MINUTES,
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
        _MINIMUM_MEMORY_MIB, _round_up(memory_request, quantum=256)
    )
    observed_output_bytes = _directory_bytes(task_record.output_root)
    estimated_output_mib = max(1, ceil(observed_output_bytes / _MEBIBYTE))
    scratch_request = 128 + ceil((observed_output_bytes * 64) / _MEBIBYTE) + ceil(
        dataset_scale * 32
    )
    estimated_scratch_mib = max(
        _MINIMUM_SCRATCH_MIB, _round_up(scratch_request, quantum=128)
    )
    estimated_core_hours = round(
        (estimated_cpus_per_task * estimated_wallclock_minutes) / 60, 2
    )
    task_log_path = Path("parallel-logs", f"{variant_run.config.variant_id}.log").as_posix()
    return RabiesMethodSensitivitySlurmJobPlanRow(
        dataset_id=dataset_id,
        variant_id=variant_run.config.variant_id,
        job_name=job_name,
        alignment_mode=variant_run.config.alignment_mode,
        trimming_mode=variant_run.config.trimming_mode,
        aligned_site_count=variant_run.alignment_length,
        trimmed_site_count=variant_run.trimmed_alignment_length,
        bootstrap_replicates=bootstrap_replicates,
        iqtree_threads=iqtree_threads,
        estimated_cpus_per_task=estimated_cpus_per_task,
        estimated_memory_mib=estimated_memory_mib,
        estimated_wallclock_minutes=estimated_wallclock_minutes,
        slurm_time=_format_slurm_time(estimated_wallclock_minutes),
        estimated_scratch_mib=estimated_scratch_mib,
        observed_output_bytes=observed_output_bytes,
        estimated_output_mib=estimated_output_mib,
        estimated_core_hours=estimated_core_hours,
        bundle_output_directory=Path("variants", variant_run.config.variant_id).as_posix(),
        task_log_path=task_log_path,
        suggested_sbatch_options=(
            f"--job-name={job_name} "
            f"--cpus-per-task={estimated_cpus_per_task} "
            f"--mem={estimated_memory_mib}M "
            f"--time={_format_slurm_time(estimated_wallclock_minutes)} "
            f"--output=slurm-logs/{variant_run.config.variant_id}.%j.out "
            f"--error=slurm-logs/{variant_run.config.variant_id}.%j.err"
        ),
    )


def _directory_bytes(path: Path) -> int:
    if not path.is_dir():
        return 0
    return sum(
        candidate.stat().st_size for candidate in path.rglob("*") if candidate.is_file()
    )


def _format_float(value: float) -> str:
    text = f"{value:.2f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def _format_slurm_time(minutes: int) -> str:
    hours, remaining_minutes = divmod(minutes, 60)
    return f"{hours:02d}:{remaining_minutes:02d}:00"


def _round_up(value: int, *, quantum: int) -> int:
    return ((value + quantum - 1) // quantum) * quantum


def _write_tsv(path: Path, rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )
    return path
