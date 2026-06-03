from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmPlanningReport
from .shared import format_float, write_tsv


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
            "taxon_count",
            "dataset_size_class",
            "method_group",
            "resource_class",
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
                str(row.taxon_count),
                row.dataset_size_class,
                row.method_group,
                row.resource_class,
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
                format_float(row.estimated_core_hours),
                row.bundle_output_directory,
                row.task_log_path,
                row.suggested_sbatch_options,
            ]
        )
    return write_tsv(path, rows)


def write_rabies_method_sensitivity_slurm_assumptions_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmPlanningReport,
) -> Path:
    """Write the sizing assumptions behind the Slurm planning report."""
    rows = [["assumption_id", "parameter", "value", "rationale"]]
    for row in report.assumptions:
        rows.append([row.assumption_id, row.parameter, row.value, row.rationale])
    return write_tsv(path, rows)


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
