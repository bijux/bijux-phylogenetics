from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm import (
    build_rabies_method_sensitivity_slurm_planning_report,
    write_rabies_method_sensitivity_slurm_assumptions_table,
    write_rabies_method_sensitivity_slurm_job_plan_table,
    write_rabies_method_sensitivity_slurm_summary_json,
)


def _write_bytes(path: Path, *, byte_count: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * byte_count)


def test_build_rabies_method_sensitivity_slurm_planning_report_sizes_jobs(
    tmp_path: Path,
) -> None:
    auto_root = tmp_path / "variants" / "auto-gap-threshold"
    ginsi_root = tmp_path / "variants" / "ginsi-gappyout"
    _write_bytes(auto_root / "alignment.aln", byte_count=48_000)
    _write_bytes(ginsi_root / "alignment.aln", byte_count=2_500_000)

    report = SimpleNamespace(
        dataset=SimpleNamespace(
            dataset_id="rabies_method_sensitivity_panel",
            workflow_prefix="rabies-method-sensitivity-panel",
            taxon_count=9,
        ),
        task_records=(
            SimpleNamespace(
                variant_id="auto-gap-threshold",
                output_root=auto_root,
            ),
            SimpleNamespace(
                variant_id="ginsi-gappyout",
                output_root=ginsi_root,
            ),
        ),
        variant_runs=(
            SimpleNamespace(
                config=SimpleNamespace(
                    variant_id="auto-gap-threshold",
                    alignment_mode="auto",
                    trimming_mode="gap-threshold",
                ),
                alignment_length=1353,
                trimmed_alignment_length=1353,
            ),
            SimpleNamespace(
                config=SimpleNamespace(
                    variant_id="ginsi-gappyout",
                    alignment_mode="ginsi",
                    trimming_mode="gappyout",
                ),
                alignment_length=1353,
                trimmed_alignment_length=1353,
            ),
        ),
        iqtree_threads=1,
        bootstrap_replicates=1000,
    )

    planning = build_rabies_method_sensitivity_slurm_planning_report(report)

    assert planning.dataset_id == "rabies_method_sensitivity_panel"
    assert planning.job_count == 2
    assert len(planning.assumptions) == 5
    auto_job, ginsi_job = planning.rows
    assert auto_job.variant_id == "auto-gap-threshold"
    assert auto_job.dataset_size_class == "compact"
    assert auto_job.method_group == "mafft-auto"
    assert auto_job.resource_class == "standard"
    assert auto_job.estimated_cpus_per_task == 1
    assert auto_job.estimated_memory_mib == 1024
    assert auto_job.estimated_wallclock_minutes >= 20
    assert auto_job.estimated_output_mib == 1
    assert "--job-name=rabies-method-sensitivity-panel-auto-gap-threshold" in (
        auto_job.suggested_sbatch_options
    )
    assert ginsi_job.variant_id == "ginsi-gappyout"
    assert ginsi_job.dataset_size_class == "compact"
    assert ginsi_job.method_group == "mafft-ginsi"
    assert ginsi_job.resource_class == "elevated"
    assert ginsi_job.estimated_cpus_per_task == 2
    assert ginsi_job.estimated_memory_mib > auto_job.estimated_memory_mib
    assert ginsi_job.estimated_wallclock_minutes > auto_job.estimated_wallclock_minutes
    assert ginsi_job.estimated_output_mib == 3
    assert ginsi_job.estimated_scratch_mib > auto_job.estimated_scratch_mib
    assert planning.total_estimated_core_hours == round(
        auto_job.estimated_core_hours + ginsi_job.estimated_core_hours,
        2,
    )
    assert planning.maximum_estimated_memory_mib == ginsi_job.estimated_memory_mib
    assert (
        planning.maximum_estimated_wallclock_minutes
        == ginsi_job.estimated_wallclock_minutes
    )


def test_write_rabies_method_sensitivity_slurm_artifacts(tmp_path: Path) -> None:
    job_root = tmp_path / "variants" / "auto-gap-threshold"
    _write_bytes(job_root / "alignment.aln", byte_count=64_000)
    report = SimpleNamespace(
        dataset=SimpleNamespace(
            dataset_id="rabies_method_sensitivity_panel",
            workflow_prefix="rabies-method-sensitivity-panel",
            taxon_count=9,
        ),
        task_records=(
            SimpleNamespace(
                variant_id="auto-gap-threshold",
                output_root=job_root,
            ),
        ),
        variant_runs=(
            SimpleNamespace(
                config=SimpleNamespace(
                    variant_id="auto-gap-threshold",
                    alignment_mode="auto",
                    trimming_mode="gap-threshold",
                ),
                alignment_length=1353,
                trimmed_alignment_length=1353,
            ),
        ),
        iqtree_threads=1,
        bootstrap_replicates=1000,
    )
    planning = build_rabies_method_sensitivity_slurm_planning_report(report)

    job_plan_path = write_rabies_method_sensitivity_slurm_job_plan_table(
        tmp_path / "slurm-job-plan.tsv",
        planning,
    )
    assumptions_path = write_rabies_method_sensitivity_slurm_assumptions_table(
        tmp_path / "slurm-estimation-assumptions.tsv",
        planning,
    )
    summary_path = write_rabies_method_sensitivity_slurm_summary_json(
        tmp_path / "slurm-planning-summary.json",
        planning,
    )

    job_plan_text = job_plan_path.read_text(encoding="utf-8")
    assert "suggested_sbatch_options" in job_plan_text
    assert "dataset_size_class" in job_plan_text
    assert "resource_class" in job_plan_text
    assert "auto-gap-threshold" in job_plan_text
    assumptions_text = assumptions_path.read_text(encoding="utf-8")
    assert "observed-output-footprint" in assumptions_text
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["dataset_id"] == "rabies_method_sensitivity_panel"
    assert payload["job_count"] == 1
    assert payload["jobs"][0]["variant_id"] == "auto-gap-threshold"
