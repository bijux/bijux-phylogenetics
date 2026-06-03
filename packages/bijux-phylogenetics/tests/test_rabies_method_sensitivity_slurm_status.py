from __future__ import annotations

from dataclasses import replace
import json
import os
from pathlib import Path

from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    load_rabies_method_sensitivity_panel_dataset,
)
from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm.status import (
    build_rabies_method_sensitivity_slurm_status_report,
    write_rabies_method_sensitivity_slurm_job_status_table,
    write_rabies_method_sensitivity_slurm_partition_status_table,
    write_rabies_method_sensitivity_slurm_status_json,
)
from bijux_phylogenetics.engines.common import (
    EngineActiveRunRecord,
    write_active_engine_run,
)


def _write_bundle_scaffold(
    root: Path,
    *,
    partition_rows: list[dict[str, str]],
    member_rows: list[dict[str, str]],
) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / "workflow-config.resolved.json").write_text(
        json.dumps(
            {
                "dataset_id": "rabies_method_sensitivity_panel",
                "workflow_prefix": "rabies-method-sensitivity-panel",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    _write_tsv(
        root / "slurm-array-partitions.tsv",
        fieldnames=(
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
        ),
        rows=partition_rows,
    )
    _write_tsv(
        root / "slurm-array-members.tsv",
        fieldnames=(
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
        ),
        rows=member_rows,
    )
    for row in partition_rows:
        script_path = root / str(row["script_path"])
        script_path.parent.mkdir(parents=True, exist_ok=True)
        script_path.write_text("#!/bin/sh\n", encoding="utf-8")


def _write_task_log(root: Path, variant_id: str, *, status: str) -> None:
    path = root / "parallel-logs" / f"{variant_id}.log"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"variant_id: {variant_id}",
                f"label: {variant_id}",
                "execution_mode: parallel",
                "alignment_mode: auto",
                "trimming_mode: gap-threshold",
                "trim_gap_threshold: 0.1",
                f"status: {status}",
                f"output_root: variants/{variant_id}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_required_variant_outputs(root: Path, variant_id: str) -> None:
    variant_root = root / "variants" / variant_id
    variant_root.mkdir(parents=True, exist_ok=True)
    for filename in (
        f"{variant_id}.aln",
        f"{variant_id}.trimmed.aln",
        "fasttree.nwk",
        "iqtree-support.nwk",
        "rooted-engine-comparison.tsv",
        "rooted-fasttree.nwk",
        "rooted-iqtree-support.nwk",
        "rooting-summary.tsv",
        "unrooted-comparison.tsv",
        "unrooted-conclusions.tsv",
        "unrooted-conflicting-clades.tsv",
        "unrooted-shared-clades.tsv",
        "unrooted-stability-summary.tsv",
        "unrooted-support-weighted-conflicts.tsv",
    ):
        (variant_root / filename).write_text("x\n", encoding="utf-8")


def _write_tsv(
    path: Path,
    *,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, str]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = "\t".join(fieldnames)
    rendered_rows = [
        "\t".join(str(row.get(field, "")) for field in fieldnames) for row in rows
    ]
    path.write_text("\n".join([header, *rendered_rows]) + "\n", encoding="utf-8")


def test_build_rabies_method_sensitivity_slurm_status_report_classifies_live_jobs(
    tmp_path: Path,
) -> None:
    root = tmp_path / "workflow"
    _write_bundle_scaffold(
        root,
        partition_rows=[
            {
                "partition_id": "compact-mafft-auto-standard",
                "dataset_size_class": "compact",
                "method_group": "mafft-auto",
                "resource_class": "standard",
                "job_count": "2",
                "array_spec": "0-1",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "variant_ids": "completed-a,failed-b",
                "trimming_modes": "gap-threshold",
                "maximum_cpus_per_task": "1",
                "maximum_memory_mib": "1024",
                "maximum_wallclock_minutes": "20",
                "total_estimated_core_hours": "0.5",
                "suggested_sbatch_command": "sbatch --array=0-1 slurm-arrays/compact-mafft-auto-standard.sbatch",
            },
            {
                "partition_id": "compact-mafft-ginsi-elevated",
                "dataset_size_class": "compact",
                "method_group": "mafft-ginsi",
                "resource_class": "elevated",
                "job_count": "1",
                "array_spec": "0-0",
                "script_path": "slurm-arrays/compact-mafft-ginsi-elevated.sbatch",
                "variant_ids": "pending-c",
                "trimming_modes": "gappyout",
                "maximum_cpus_per_task": "2",
                "maximum_memory_mib": "1536",
                "maximum_wallclock_minutes": "35",
                "total_estimated_core_hours": "0.7",
                "suggested_sbatch_command": "sbatch --array=0-0 slurm-arrays/compact-mafft-ginsi-elevated.sbatch",
            },
        ],
        member_rows=[
            {
                "partition_id": "compact-mafft-auto-standard",
                "array_index": "0",
                "variant_id": "completed-a",
                "dataset_size_class": "compact",
                "method_group": "mafft-auto",
                "trimming_mode": "gap-threshold",
                "resource_class": "standard",
                "estimated_cpus_per_task": "1",
                "estimated_memory_mib": "1024",
                "estimated_wallclock_minutes": "20",
                "estimated_core_hours": "0.2",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "bundle_output_directory": "variants/completed-a",
                "task_log_path": "parallel-logs/completed-a.log",
            },
            {
                "partition_id": "compact-mafft-auto-standard",
                "array_index": "1",
                "variant_id": "failed-b",
                "dataset_size_class": "compact",
                "method_group": "mafft-auto",
                "trimming_mode": "gap-threshold",
                "resource_class": "standard",
                "estimated_cpus_per_task": "1",
                "estimated_memory_mib": "1024",
                "estimated_wallclock_minutes": "20",
                "estimated_core_hours": "0.3",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "bundle_output_directory": "variants/failed-b",
                "task_log_path": "parallel-logs/failed-b.log",
            },
            {
                "partition_id": "compact-mafft-ginsi-elevated",
                "array_index": "0",
                "variant_id": "pending-c",
                "dataset_size_class": "compact",
                "method_group": "mafft-ginsi",
                "trimming_mode": "gappyout",
                "resource_class": "elevated",
                "estimated_cpus_per_task": "2",
                "estimated_memory_mib": "1536",
                "estimated_wallclock_minutes": "35",
                "estimated_core_hours": "0.7",
                "script_path": "slurm-arrays/compact-mafft-ginsi-elevated.sbatch",
                "bundle_output_directory": "variants/pending-c",
                "task_log_path": "parallel-logs/pending-c.log",
            },
        ],
    )
    _write_task_log(root, "completed-a", status="succeeded")
    _write_task_log(root, "failed-b", status="failed")
    _write_required_variant_outputs(root, "completed-a")
    write_active_engine_run(
        EngineActiveRunRecord(
            engine_name="bijux",
            workflow="rabies_method_sensitivity_panel",
            executable="bijux-phylogenetics",
            working_directory=root,
            manifest_path=root / "rabies-method-sensitivity-panel.run.json",
            command=["run_rabies_method_sensitivity_panel_workflow"],
            process_id=os.getpid(),
            started_at_utc="2026-05-18T12:00:00Z",
        )
    )

    report = build_rabies_method_sensitivity_slurm_status_report(root)

    assert report.active_run_state == "live"
    assert report.job_count == 3
    assert report.completed_job_count == 1
    assert report.failed_job_count == 1
    assert report.pending_job_count == 1
    assert report.stale_job_count == 0
    assert [row.status for row in report.jobs] == ["completed", "failed", "pending"]
    assert report.partitions[0].overall_status == "mixed"
    assert report.partitions[1].overall_status == "pending"


def test_build_rabies_method_sensitivity_slurm_status_report_marks_stale_dead_marker(
    tmp_path: Path,
) -> None:
    root = tmp_path / "workflow"
    _write_bundle_scaffold(
        root,
        partition_rows=[
            {
                "partition_id": "compact-mafft-auto-standard",
                "dataset_size_class": "compact",
                "method_group": "mafft-auto",
                "resource_class": "standard",
                "job_count": "1",
                "array_spec": "0-0",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "variant_ids": "stale-a",
                "trimming_modes": "gap-threshold",
                "maximum_cpus_per_task": "1",
                "maximum_memory_mib": "1024",
                "maximum_wallclock_minutes": "20",
                "total_estimated_core_hours": "0.2",
                "suggested_sbatch_command": "sbatch --array=0-0 slurm-arrays/compact-mafft-auto-standard.sbatch",
            }
        ],
        member_rows=[
            {
                "partition_id": "compact-mafft-auto-standard",
                "array_index": "0",
                "variant_id": "stale-a",
                "dataset_size_class": "compact",
                "method_group": "mafft-auto",
                "trimming_mode": "gap-threshold",
                "resource_class": "standard",
                "estimated_cpus_per_task": "1",
                "estimated_memory_mib": "1024",
                "estimated_wallclock_minutes": "20",
                "estimated_core_hours": "0.2",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "bundle_output_directory": "variants/stale-a",
                "task_log_path": "parallel-logs/stale-a.log",
            }
        ],
    )
    write_active_engine_run(
        EngineActiveRunRecord(
            engine_name="bijux",
            workflow="rabies_method_sensitivity_panel",
            executable="bijux-phylogenetics",
            working_directory=root,
            manifest_path=root / "rabies-method-sensitivity-panel.run.json",
            command=["run_rabies_method_sensitivity_panel_workflow"],
            process_id=0,
            started_at_utc="2026-05-18T12:00:00Z",
        )
    )

    report = build_rabies_method_sensitivity_slurm_status_report(root)

    assert report.active_run_state == "stale"
    assert report.jobs[0].status == "stale"
    assert report.jobs[0].evidence_class == "stale-running-marker"
    assert report.stale_job_count == 1


def test_build_rabies_method_sensitivity_slurm_status_report_marks_partial_outputs_stale(
    tmp_path: Path,
) -> None:
    root = tmp_path / "workflow"
    _write_bundle_scaffold(
        root,
        partition_rows=[
            {
                "partition_id": "compact-mafft-auto-standard",
                "dataset_size_class": "compact",
                "method_group": "mafft-auto",
                "resource_class": "standard",
                "job_count": "1",
                "array_spec": "0-0",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "variant_ids": "stale-b",
                "trimming_modes": "gap-threshold",
                "maximum_cpus_per_task": "1",
                "maximum_memory_mib": "1024",
                "maximum_wallclock_minutes": "20",
                "total_estimated_core_hours": "0.2",
                "suggested_sbatch_command": "sbatch --array=0-0 slurm-arrays/compact-mafft-auto-standard.sbatch",
            }
        ],
        member_rows=[
            {
                "partition_id": "compact-mafft-auto-standard",
                "array_index": "0",
                "variant_id": "stale-b",
                "dataset_size_class": "compact",
                "method_group": "mafft-auto",
                "trimming_mode": "gap-threshold",
                "resource_class": "standard",
                "estimated_cpus_per_task": "1",
                "estimated_memory_mib": "1024",
                "estimated_wallclock_minutes": "20",
                "estimated_core_hours": "0.2",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "bundle_output_directory": "variants/stale-b",
                "task_log_path": "parallel-logs/stale-b.log",
            }
        ],
    )
    variant_root = root / "variants" / "stale-b"
    variant_root.mkdir(parents=True, exist_ok=True)
    (variant_root / "partial.txt").write_text("partial\n", encoding="utf-8")

    report = build_rabies_method_sensitivity_slurm_status_report(root)

    assert report.active_run_state == "absent"
    assert report.jobs[0].status == "stale"
    assert report.jobs[0].evidence_class == "abandoned-partial-output"
    assert report.jobs[0].output_file_count == 1
    assert report.jobs[0].missing_required_file_count > 0


def test_build_rabies_method_sensitivity_slurm_status_report_marks_completed_jobs_stale_when_settings_drift() -> (
    None
):
    dataset = load_rabies_method_sensitivity_panel_dataset()
    changed_variant = replace(
        dataset.variants[0],
        trim_gap_threshold=dataset.variants[0].trim_gap_threshold + 0.05,
    )
    drifted_dataset = replace(
        dataset,
        variants=(changed_variant, *dataset.variants[1:]),
    )

    report = build_rabies_method_sensitivity_slurm_status_report(
        dataset.reference_output_root,
        dataset=drifted_dataset,
    )

    stale_rows = [row for row in report.jobs if row.status == "stale"]
    assert [row.variant_id for row in stale_rows] == [dataset.variants[0].variant_id]
    assert stale_rows[0].evidence_class == "stale-output-drift"
    assert stale_rows[0].output_freshness_status == "stale"
    assert report.stale_job_count == 1
    assert report.stale_output_job_count == 1


def test_write_rabies_method_sensitivity_slurm_status_artifacts(tmp_path: Path) -> None:
    root = tmp_path / "workflow"
    _write_bundle_scaffold(
        root,
        partition_rows=[
            {
                "partition_id": "compact-mafft-auto-standard",
                "dataset_size_class": "compact",
                "method_group": "mafft-auto",
                "resource_class": "standard",
                "job_count": "1",
                "array_spec": "0-0",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "variant_ids": "completed-a",
                "trimming_modes": "gap-threshold",
                "maximum_cpus_per_task": "1",
                "maximum_memory_mib": "1024",
                "maximum_wallclock_minutes": "20",
                "total_estimated_core_hours": "0.2",
                "suggested_sbatch_command": "sbatch --array=0-0 slurm-arrays/compact-mafft-auto-standard.sbatch",
            }
        ],
        member_rows=[
            {
                "partition_id": "compact-mafft-auto-standard",
                "array_index": "0",
                "variant_id": "completed-a",
                "dataset_size_class": "compact",
                "method_group": "mafft-auto",
                "trimming_mode": "gap-threshold",
                "resource_class": "standard",
                "estimated_cpus_per_task": "1",
                "estimated_memory_mib": "1024",
                "estimated_wallclock_minutes": "20",
                "estimated_core_hours": "0.2",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "bundle_output_directory": "variants/completed-a",
                "task_log_path": "parallel-logs/completed-a.log",
            }
        ],
    )
    _write_task_log(root, "completed-a", status="succeeded")
    _write_required_variant_outputs(root, "completed-a")

    report = build_rabies_method_sensitivity_slurm_status_report(root)
    job_status_path = write_rabies_method_sensitivity_slurm_job_status_table(
        tmp_path / "slurm-job-status.tsv",
        report,
    )
    partition_status_path = (
        write_rabies_method_sensitivity_slurm_partition_status_table(
            tmp_path / "slurm-partition-status.tsv",
            report,
        )
    )
    workflow_status_path = write_rabies_method_sensitivity_slurm_status_json(
        tmp_path / "slurm-workflow-status.json",
        report,
    )

    assert "missing_required_file_count" in job_status_path.read_text(encoding="utf-8")
    assert "output_freshness_status" in job_status_path.read_text(encoding="utf-8")
    assert "overall_status" in partition_status_path.read_text(encoding="utf-8")
    payload = json.loads(workflow_status_path.read_text(encoding="utf-8"))
    assert payload["bundle_root"] == "."
    assert (
        payload["execution_record_path"] == "rabies-method-sensitivity-panel.run.json"
    )
    assert payload["completed_job_count"] == 1
    assert payload["fresh_output_job_count"] == 1
    assert payload["jobs"][0]["status"] == "completed"
