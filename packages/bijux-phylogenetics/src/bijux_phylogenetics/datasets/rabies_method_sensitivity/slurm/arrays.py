from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path

from .planning import (
    RabiesMethodSensitivitySlurmJobPlanRow,
    RabiesMethodSensitivitySlurmPlanningReport,
)

__all__ = [
    "RabiesMethodSensitivitySlurmArrayMemberRow",
    "RabiesMethodSensitivitySlurmArrayPartitionRow",
    "RabiesMethodSensitivitySlurmArrayStrategyReport",
    "build_rabies_method_sensitivity_slurm_array_strategy_report",
    "write_rabies_method_sensitivity_slurm_array_members_table",
    "write_rabies_method_sensitivity_slurm_array_partition_scripts",
    "write_rabies_method_sensitivity_slurm_array_partitions_table",
    "write_rabies_method_sensitivity_slurm_array_strategy_json",
]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmArrayMemberRow:
    """One job assigned to one concrete Slurm array partition."""

    partition_id: str
    array_index: int
    variant_id: str
    dataset_size_class: str
    method_group: str
    trimming_mode: str
    resource_class: str
    estimated_cpus_per_task: int
    estimated_memory_mib: int
    estimated_wallclock_minutes: int
    estimated_core_hours: float
    script_path: str
    bundle_output_directory: str
    task_log_path: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmArrayPartitionRow:
    """One array-ready partition over jobs that can share one envelope."""

    partition_id: str
    dataset_size_class: str
    method_group: str
    resource_class: str
    job_count: int
    array_spec: str
    script_path: str
    variant_ids: tuple[str, ...]
    trimming_modes: tuple[str, ...]
    maximum_cpus_per_task: int
    maximum_memory_mib: int
    maximum_wallclock_minutes: int
    total_estimated_core_hours: float
    suggested_sbatch_command: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmArrayStrategyReport:
    """One reviewer-facing partition strategy for Slurm job arrays."""

    dataset_id: str
    workflow_prefix: str
    partition_count: int
    total_job_count: int
    script_count: int
    largest_partition_size: int
    partitions: tuple[RabiesMethodSensitivitySlurmArrayPartitionRow, ...]
    members: tuple[RabiesMethodSensitivitySlurmArrayMemberRow, ...]


def build_rabies_method_sensitivity_slurm_array_strategy_report(
    planning_report: RabiesMethodSensitivitySlurmPlanningReport,
) -> RabiesMethodSensitivitySlurmArrayStrategyReport:
    """Group Slurm jobs into real array partitions by size, method, and resources."""
    grouped_jobs: dict[
        tuple[str, str, str],
        list[RabiesMethodSensitivitySlurmJobPlanRow],
    ] = {}
    for row in planning_report.rows:
        key = (row.dataset_size_class, row.method_group, row.resource_class)
        grouped_jobs.setdefault(key, []).append(row)

    partition_rows: list[RabiesMethodSensitivitySlurmArrayPartitionRow] = []
    member_rows: list[RabiesMethodSensitivitySlurmArrayMemberRow] = []
    for dataset_size_class, method_group, resource_class in sorted(grouped_jobs):
        partition_id = f"{dataset_size_class}-{method_group}-{resource_class}"
        jobs = sorted(
            grouped_jobs[(dataset_size_class, method_group, resource_class)],
            key=lambda row: row.variant_id,
        )
        array_spec = f"0-{len(jobs) - 1}"
        script_path = Path("slurm-arrays", f"{partition_id}.sbatch").as_posix()
        partition_rows.append(
            RabiesMethodSensitivitySlurmArrayPartitionRow(
                partition_id=partition_id,
                dataset_size_class=dataset_size_class,
                method_group=method_group,
                resource_class=resource_class,
                job_count=len(jobs),
                array_spec=array_spec,
                script_path=script_path,
                variant_ids=tuple(job.variant_id for job in jobs),
                trimming_modes=tuple(sorted({job.trimming_mode for job in jobs})),
                maximum_cpus_per_task=max(job.estimated_cpus_per_task for job in jobs),
                maximum_memory_mib=max(job.estimated_memory_mib for job in jobs),
                maximum_wallclock_minutes=max(
                    job.estimated_wallclock_minutes for job in jobs
                ),
                total_estimated_core_hours=round(
                    sum(job.estimated_core_hours for job in jobs),
                    2,
                ),
                suggested_sbatch_command=(f"sbatch --array={array_spec} {script_path}"),
            )
        )
        for array_index, job in enumerate(jobs):
            member_rows.append(
                RabiesMethodSensitivitySlurmArrayMemberRow(
                    partition_id=partition_id,
                    array_index=array_index,
                    variant_id=job.variant_id,
                    dataset_size_class=dataset_size_class,
                    method_group=method_group,
                    trimming_mode=job.trimming_mode,
                    resource_class=resource_class,
                    estimated_cpus_per_task=job.estimated_cpus_per_task,
                    estimated_memory_mib=job.estimated_memory_mib,
                    estimated_wallclock_minutes=job.estimated_wallclock_minutes,
                    estimated_core_hours=job.estimated_core_hours,
                    script_path=script_path,
                    bundle_output_directory=job.bundle_output_directory,
                    task_log_path=job.task_log_path,
                )
            )

    return RabiesMethodSensitivitySlurmArrayStrategyReport(
        dataset_id=planning_report.dataset_id,
        workflow_prefix=planning_report.workflow_prefix,
        partition_count=len(partition_rows),
        total_job_count=len(planning_report.rows),
        script_count=len(partition_rows),
        largest_partition_size=max(row.job_count for row in partition_rows),
        partitions=tuple(partition_rows),
        members=tuple(member_rows),
    )


def write_rabies_method_sensitivity_slurm_array_partitions_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmArrayStrategyReport,
) -> Path:
    """Write one partition-level array strategy ledger."""
    rows = [
        [
            "partition_id",
            "dataset_size_class",
            "method_group",
            "resource_class",
            "job_count",
            "array_spec",
            "script_path",
            "variant_ids",
            "trimming_modes",
            "maximum_cpus_per_task",
            "maximum_memory_mib",
            "maximum_wallclock_minutes",
            "total_estimated_core_hours",
            "suggested_sbatch_command",
        ]
    ]
    for row in report.partitions:
        rows.append(
            [
                row.partition_id,
                row.dataset_size_class,
                row.method_group,
                row.resource_class,
                str(row.job_count),
                row.array_spec,
                row.script_path,
                ",".join(row.variant_ids),
                ",".join(row.trimming_modes),
                str(row.maximum_cpus_per_task),
                str(row.maximum_memory_mib),
                str(row.maximum_wallclock_minutes),
                _format_float(row.total_estimated_core_hours),
                row.suggested_sbatch_command,
            ]
        )
    return _write_tsv(path, rows)


def write_rabies_method_sensitivity_slurm_array_members_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmArrayStrategyReport,
) -> Path:
    """Write one member-level Slurm array assignment ledger."""
    rows = [
        [
            "partition_id",
            "array_index",
            "variant_id",
            "dataset_size_class",
            "method_group",
            "trimming_mode",
            "resource_class",
            "estimated_cpus_per_task",
            "estimated_memory_mib",
            "estimated_wallclock_minutes",
            "estimated_core_hours",
            "script_path",
            "bundle_output_directory",
            "task_log_path",
        ]
    ]
    for row in report.members:
        rows.append(
            [
                row.partition_id,
                str(row.array_index),
                row.variant_id,
                row.dataset_size_class,
                row.method_group,
                row.trimming_mode,
                row.resource_class,
                str(row.estimated_cpus_per_task),
                str(row.estimated_memory_mib),
                str(row.estimated_wallclock_minutes),
                _format_float(row.estimated_core_hours),
                row.script_path,
                row.bundle_output_directory,
                row.task_log_path,
            ]
        )
    return _write_tsv(path, rows)


def write_rabies_method_sensitivity_slurm_array_strategy_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmArrayStrategyReport,
) -> Path:
    """Write the structured array partition strategy."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "dataset_id": report.dataset_id,
        "workflow_prefix": report.workflow_prefix,
        "partition_count": report.partition_count,
        "total_job_count": report.total_job_count,
        "script_count": report.script_count,
        "largest_partition_size": report.largest_partition_size,
        "partitions": [asdict(row) for row in report.partitions],
        "members": [asdict(row) for row in report.members],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_rabies_method_sensitivity_slurm_array_partition_scripts(
    output_root: Path,
    report: RabiesMethodSensitivitySlurmArrayStrategyReport,
) -> Path:
    """Write one executable sbatch script per partition."""
    output_root.mkdir(parents=True, exist_ok=True)
    members_by_partition: dict[
        str, list[RabiesMethodSensitivitySlurmArrayMemberRow]
    ] = {}
    for member in report.members:
        members_by_partition.setdefault(member.partition_id, []).append(member)
    for partition in report.partitions:
        members = sorted(
            members_by_partition[partition.partition_id],
            key=lambda member: member.array_index,
        )
        script_path = output_root / Path(partition.script_path).name
        script_lines = [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "",
            f"# partition_id: {partition.partition_id}",
            f"# dataset_size_class: {partition.dataset_size_class}",
            f"# method_group: {partition.method_group}",
            f"# resource_class: {partition.resource_class}",
            f"# array_spec: {partition.array_spec}",
            f"# maximum_cpus_per_task: {partition.maximum_cpus_per_task}",
            f"# maximum_memory_mib: {partition.maximum_memory_mib}",
            f"# maximum_wallclock_minutes: {partition.maximum_wallclock_minutes}",
            "",
            f"#SBATCH --job-name={partition.partition_id}",
            f"#SBATCH --cpus-per-task={partition.maximum_cpus_per_task}",
            f"#SBATCH --mem={partition.maximum_memory_mib}M",
            (
                "#SBATCH --time="
                f"{_format_slurm_time(partition.maximum_wallclock_minutes)}"
            ),
            f"#SBATCH --array={partition.array_spec}",
            "",
            "variant_ids=(",
        ]
        script_lines.extend(f'  "{member.variant_id}"' for member in members)
        script_lines.extend(
            [
                ")",
                "",
                'if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then',
                '  echo "SLURM_ARRAY_TASK_ID is required" >&2',
                "  exit 2",
                "fi",
                "",
                'variant_id="${variant_ids[$SLURM_ARRAY_TASK_ID]}"',
                (
                    'array_root="${BIJUX_PHYLOGENETICS_ARRAY_ROOT:-artifacts/'
                    'rabies-method-sensitivity-arrays}"'
                ),
                'output_root="${array_root}/'
                f'{partition.partition_id}/${{variant_id}}"',
                "exec bijux-phylogenetics demo rabies-method-sensitivity-panel \\",
                '  --out "${output_root}" \\',
                '  --variant-id "${variant_id}" \\',
                "  --parallel-workers 1 \\",
                '  "$@"',
            ]
        )
        script_path.write_text("\n".join(script_lines) + "\n", encoding="utf-8")
        script_path.chmod(0o755)
    return output_root


def _format_float(value: float) -> str:
    text = f"{value:.2f}"
    return text.rstrip("0").rstrip(".") if "." in text else text


def _format_slurm_time(minutes: int) -> str:
    hours, remaining_minutes = divmod(minutes, 60)
    return f"{hours:02d}:{remaining_minutes:02d}:00"


def _write_tsv(path: Path, rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )
    return path
