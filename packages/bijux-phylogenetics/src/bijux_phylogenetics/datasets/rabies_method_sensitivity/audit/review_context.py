from __future__ import annotations

from dataclasses import dataclass

from .inventory import RabiesMethodSensitivityAuditSnapshot
from .io import _parse_task_log


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityAuditReviewContext:
    config_variants: dict[str, dict[str, object]]
    manifest_task_records: dict[str, dict[str, object]]
    parallel_summary_rows: dict[str, dict[str, str]]
    variant_summary_rows: dict[str, dict[str, str]]
    logged_variants: dict[str, dict[str, str]]
    config_variant_ids: list[str]
    manifest_variant_ids: list[str]
    parallel_variant_ids: list[str]
    summary_variant_ids: list[str]
    logged_variant_ids: list[str]
    written_variant_ids: list[str]
    expected_variant_count: int
    slurm_partition_ids: set[str]
    slurm_script_paths: set[str]
    member_partition_ids: set[str]
    member_variant_ids: list[str]
    slurm_job_evidence_variant_ids: list[str]
    slurm_storage_category_ids: list[str]
    slurm_storage_variant_ids: list[str]
    slurm_output_explosion_variant_ids: list[str]
    slurm_tree_retention_variant_ids: list[str]
    slurm_merge_variant_ids: list[str]
    slurm_output_freshness_variant_ids: list[str]
    slurm_job_status_variant_ids: list[str]
    slurm_partition_status_ids: list[str]
    slurm_output_freshness_rows_by_variant: dict[str, dict[str, str]]
    slurm_job_evidence_rows_by_variant: dict[str, dict[str, str]]
    slurm_merge_rows_by_variant: dict[str, dict[str, str]]
    slurm_job_status_rows_by_variant: dict[str, dict[str, str]]
    slurm_failure_recovery_variant_ids: list[str]


def build_rabies_method_sensitivity_audit_review_context(
    snapshot: RabiesMethodSensitivityAuditSnapshot,
) -> RabiesMethodSensitivityAuditReviewContext:
    config_variants = {
        str(row["variant_id"]): row
        for row in list(snapshot.resolved_config.get("variants", []))
    }
    manifest_task_records = {
        str(row["variant_id"]): row
        for row in list(snapshot.workflow_manifest.get("task_records", []))
    }
    parallel_summary_rows = {
        str(row["variant_id"]): row for row in snapshot.parallel_rows
    }
    variant_summary_rows = {
        str(row["variant_id"]): row for row in snapshot.variant_rows
    }
    logged_variants = {
        path.stem: _parse_task_log(path)
        for path in sorted(snapshot.task_logs_root.glob("*.log"))
    }
    config_variant_ids = sorted(config_variants)
    manifest_variant_ids = sorted(manifest_task_records)
    parallel_variant_ids = sorted(parallel_summary_rows)
    summary_variant_ids = sorted(variant_summary_rows)
    logged_variant_ids = sorted(logged_variants)
    written_variant_ids = sorted(
        path.name for path in snapshot.variants_root.iterdir() if path.is_dir()
    )
    slurm_partition_ids = {
        str(row["partition_id"]) for row in snapshot.slurm_array_partition_rows
    }
    slurm_script_paths = {
        str(row["script_path"]) for row in snapshot.slurm_array_partition_rows
    }
    member_partition_ids = {
        str(row["partition_id"]) for row in snapshot.slurm_array_member_rows
    }
    member_variant_ids = sorted(
        str(row["variant_id"]) for row in snapshot.slurm_array_member_rows
    )
    slurm_job_evidence_variant_ids = sorted(
        str(row["variant_id"]) for row in snapshot.slurm_job_evidence_rows
    )
    slurm_storage_category_ids = sorted(
        str(row["category_id"]) for row in snapshot.slurm_storage_category_rows
    )
    slurm_storage_variant_ids = sorted(
        str(row["variant_id"]) for row in snapshot.slurm_storage_variant_rows
    )
    slurm_output_explosion_variant_ids = sorted(
        str(row["variant_id"]) for row in snapshot.slurm_output_explosion_variant_rows
    )
    slurm_tree_retention_variant_ids = sorted(
        {str(row["variant_id"]) for row in snapshot.slurm_tree_retention_file_rows}
    )
    slurm_merge_variant_ids = sorted(
        str(row["variant_id"]) for row in snapshot.slurm_merge_variant_rows
    )
    slurm_output_freshness_variant_ids = sorted(
        str(row["variant_id"]) for row in snapshot.slurm_output_freshness_rows
    )
    slurm_job_status_variant_ids = sorted(
        str(row["variant_id"]) for row in snapshot.slurm_job_status_rows
    )
    slurm_partition_status_ids = sorted(
        str(row["partition_id"]) for row in snapshot.slurm_partition_status_rows
    )
    slurm_output_freshness_rows_by_variant = {
        str(row["variant_id"]): row for row in snapshot.slurm_output_freshness_rows
    }
    slurm_job_evidence_rows_by_variant = {
        str(row["variant_id"]): row for row in snapshot.slurm_job_evidence_rows
    }
    slurm_merge_rows_by_variant = {
        str(row["variant_id"]): row for row in snapshot.slurm_merge_variant_rows
    }
    slurm_job_status_rows_by_variant = {
        str(row["variant_id"]): row for row in snapshot.slurm_job_status_rows
    }
    slurm_failure_recovery_variant_ids = sorted(
        {str(row["variant_id"]) for row in snapshot.slurm_failure_recovery_job_rows}
    )
    return RabiesMethodSensitivityAuditReviewContext(
        config_variants=config_variants,
        manifest_task_records=manifest_task_records,
        parallel_summary_rows=parallel_summary_rows,
        variant_summary_rows=variant_summary_rows,
        logged_variants=logged_variants,
        config_variant_ids=config_variant_ids,
        manifest_variant_ids=manifest_variant_ids,
        parallel_variant_ids=parallel_variant_ids,
        summary_variant_ids=summary_variant_ids,
        logged_variant_ids=logged_variant_ids,
        written_variant_ids=written_variant_ids,
        expected_variant_count=len(config_variant_ids),
        slurm_partition_ids=slurm_partition_ids,
        slurm_script_paths=slurm_script_paths,
        member_partition_ids=member_partition_ids,
        member_variant_ids=member_variant_ids,
        slurm_job_evidence_variant_ids=slurm_job_evidence_variant_ids,
        slurm_storage_category_ids=slurm_storage_category_ids,
        slurm_storage_variant_ids=slurm_storage_variant_ids,
        slurm_output_explosion_variant_ids=slurm_output_explosion_variant_ids,
        slurm_tree_retention_variant_ids=slurm_tree_retention_variant_ids,
        slurm_merge_variant_ids=slurm_merge_variant_ids,
        slurm_output_freshness_variant_ids=slurm_output_freshness_variant_ids,
        slurm_job_status_variant_ids=slurm_job_status_variant_ids,
        slurm_partition_status_ids=slurm_partition_status_ids,
        slurm_output_freshness_rows_by_variant=slurm_output_freshness_rows_by_variant,
        slurm_job_evidence_rows_by_variant=slurm_job_evidence_rows_by_variant,
        slurm_merge_rows_by_variant=slurm_merge_rows_by_variant,
        slurm_job_status_rows_by_variant=slurm_job_status_rows_by_variant,
        slurm_failure_recovery_variant_ids=slurm_failure_recovery_variant_ids,
    )
