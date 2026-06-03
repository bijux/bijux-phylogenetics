from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .io import _load_json, _read_tsv_rows

_WORKFLOW_MANIFEST_FILENAME = "rabies-method-sensitivity.manifest.json"
_REPORT_MANIFEST_FILENAME = (
    "report-artifacts/rabies-method-sensitivity-report.manifest.json"
)
_CONFIG_FILENAME = "workflow-config.resolved.json"
_PARALLEL_SUMMARY_FILENAME = "parallel-execution-summary.tsv"
_VARIANT_SUMMARY_FILENAME = "variant-summary.tsv"
_SLURM_ARRAY_PARTITIONS_FILENAME = "slurm-array-partitions.tsv"
_SLURM_ARRAY_MEMBERS_FILENAME = "slurm-array-members.tsv"
_SLURM_JOB_EVIDENCE_DIRECTORY = "slurm-job-evidence"
_SLURM_JOB_EVIDENCE_INDEX_FILENAME = "slurm-job-evidence.tsv"
_SLURM_JOB_EVIDENCE_SUMMARY_FILENAME = "slurm-job-evidence-summary.json"
_SLURM_STORAGE_CATEGORIES_FILENAME = "slurm-storage-categories.tsv"
_SLURM_STORAGE_VARIANTS_FILENAME = "slurm-storage-variants.tsv"
_SLURM_STORAGE_SUMMARY_FILENAME = "slurm-storage-report.json"
_SLURM_OUTPUT_EXPLOSION_CHECKS_FILENAME = "slurm-output-explosion-checks.tsv"
_SLURM_OUTPUT_EXPLOSION_VARIANTS_FILENAME = "slurm-output-explosion-variants.tsv"
_SLURM_OUTPUT_EXPLOSION_SUMMARY_FILENAME = "slurm-output-explosion-report.json"
_SLURM_TREE_RETENTION_CHECKS_FILENAME = "slurm-tree-retention-checks.tsv"
_SLURM_TREE_RETENTION_FILES_FILENAME = "slurm-tree-retention-files.tsv"
_SLURM_TREE_RETENTION_SUMMARY_FILENAME = "slurm-tree-retention-policy.json"
_SLURM_MERGE_CHECKS_FILENAME = "slurm-merge-checks.tsv"
_SLURM_MERGE_VARIANTS_FILENAME = "slurm-merge-variants.tsv"
_SLURM_MERGE_SUMMARY_FILENAME = "slurm-merge-report.json"
_SLURM_OUTPUT_FRESHNESS_FILENAME = "slurm-output-freshness.tsv"
_SLURM_OUTPUT_FRESHNESS_CHECKS_FILENAME = "slurm-output-freshness-checks.tsv"
_SLURM_OUTPUT_FRESHNESS_SUMMARY_FILENAME = "slurm-output-freshness.json"
_SLURM_JOB_STATUS_FILENAME = "slurm-job-status.tsv"
_SLURM_PARTITION_STATUS_FILENAME = "slurm-partition-status.tsv"
_SLURM_WORKFLOW_STATUS_FILENAME = "slurm-workflow-status.json"
_SLURM_FAILURE_RECOVERY_JOBS_FILENAME = "slurm-failure-recovery-jobs.tsv"
_SLURM_FAILURE_RECOVERY_PARTITIONS_FILENAME = "slurm-failure-recovery-partitions.tsv"
_SLURM_FAILURE_RECOVERY_SUMMARY_FILENAME = "slurm-failure-recovery-report.json"
_TASK_LOGS_DIRECTORY = "parallel-logs"
_VARIANTS_DIRECTORY = "variants"
_EXPECTED_VARIANT_FILENAMES = (
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
)


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityAuditSnapshot:
    bundle_root: Path
    workflow_manifest_path: Path
    report_manifest_path: Path
    config_path: Path
    task_logs_root: Path
    variants_root: Path
    slurm_job_evidence_root: Path
    workflow_manifest: dict[str, object]
    report_manifest: dict[str, object]
    resolved_config: dict[str, object]
    parallel_rows: list[dict[str, str]]
    variant_rows: list[dict[str, str]]
    slurm_array_partition_rows: list[dict[str, str]]
    slurm_array_member_rows: list[dict[str, str]]
    slurm_job_evidence_rows: list[dict[str, str]]
    slurm_job_evidence_summary: dict[str, object]
    slurm_storage_category_rows: list[dict[str, str]]
    slurm_storage_variant_rows: list[dict[str, str]]
    slurm_storage_summary: dict[str, object]
    slurm_output_explosion_check_rows: list[dict[str, str]]
    slurm_output_explosion_variant_rows: list[dict[str, str]]
    slurm_output_explosion_summary: dict[str, object]
    slurm_tree_retention_check_rows: list[dict[str, str]]
    slurm_tree_retention_file_rows: list[dict[str, str]]
    slurm_tree_retention_summary: dict[str, object]
    slurm_merge_check_rows: list[dict[str, str]]
    slurm_merge_variant_rows: list[dict[str, str]]
    slurm_merge_summary: dict[str, object]
    slurm_output_freshness_rows: list[dict[str, str]]
    slurm_output_freshness_check_rows: list[dict[str, str]]
    slurm_output_freshness_summary: dict[str, object]
    slurm_job_status_rows: list[dict[str, str]]
    slurm_partition_status_rows: list[dict[str, str]]
    slurm_workflow_status: dict[str, object]
    slurm_failure_recovery_job_rows: list[dict[str, str]]
    slurm_failure_recovery_partition_rows: list[dict[str, str]]
    slurm_failure_recovery_summary: dict[str, object]


def load_rabies_method_sensitivity_audit_snapshot(
    bundle_root: Path,
) -> RabiesMethodSensitivityAuditSnapshot:
    bundle_root = bundle_root.resolve()
    workflow_manifest_path = bundle_root / _WORKFLOW_MANIFEST_FILENAME
    report_manifest_path = bundle_root / _REPORT_MANIFEST_FILENAME
    config_path = bundle_root / _CONFIG_FILENAME
    task_logs_root = bundle_root / _TASK_LOGS_DIRECTORY
    variants_root = bundle_root / _VARIANTS_DIRECTORY
    slurm_job_evidence_root = bundle_root / _SLURM_JOB_EVIDENCE_DIRECTORY

    return RabiesMethodSensitivityAuditSnapshot(
        bundle_root=bundle_root,
        workflow_manifest_path=workflow_manifest_path,
        report_manifest_path=report_manifest_path,
        config_path=config_path,
        task_logs_root=task_logs_root,
        variants_root=variants_root,
        slurm_job_evidence_root=slurm_job_evidence_root,
        workflow_manifest=_load_json(workflow_manifest_path),
        report_manifest=_load_json(report_manifest_path),
        resolved_config=_load_json(config_path),
        parallel_rows=_read_tsv_rows(bundle_root / _PARALLEL_SUMMARY_FILENAME),
        variant_rows=_read_tsv_rows(bundle_root / _VARIANT_SUMMARY_FILENAME),
        slurm_array_partition_rows=_read_tsv_rows(
            bundle_root / _SLURM_ARRAY_PARTITIONS_FILENAME
        ),
        slurm_array_member_rows=_read_tsv_rows(
            bundle_root / _SLURM_ARRAY_MEMBERS_FILENAME
        ),
        slurm_job_evidence_rows=_read_tsv_rows(
            bundle_root / _SLURM_JOB_EVIDENCE_INDEX_FILENAME
        ),
        slurm_job_evidence_summary=_load_json(
            bundle_root / _SLURM_JOB_EVIDENCE_SUMMARY_FILENAME
        ),
        slurm_storage_category_rows=_read_tsv_rows(
            bundle_root / _SLURM_STORAGE_CATEGORIES_FILENAME
        ),
        slurm_storage_variant_rows=_read_tsv_rows(
            bundle_root / _SLURM_STORAGE_VARIANTS_FILENAME
        ),
        slurm_storage_summary=_load_json(bundle_root / _SLURM_STORAGE_SUMMARY_FILENAME),
        slurm_output_explosion_check_rows=_read_tsv_rows(
            bundle_root / _SLURM_OUTPUT_EXPLOSION_CHECKS_FILENAME
        ),
        slurm_output_explosion_variant_rows=_read_tsv_rows(
            bundle_root / _SLURM_OUTPUT_EXPLOSION_VARIANTS_FILENAME
        ),
        slurm_output_explosion_summary=_load_json(
            bundle_root / _SLURM_OUTPUT_EXPLOSION_SUMMARY_FILENAME
        ),
        slurm_tree_retention_check_rows=_read_tsv_rows(
            bundle_root / _SLURM_TREE_RETENTION_CHECKS_FILENAME
        ),
        slurm_tree_retention_file_rows=_read_tsv_rows(
            bundle_root / _SLURM_TREE_RETENTION_FILES_FILENAME
        ),
        slurm_tree_retention_summary=_load_json(
            bundle_root / _SLURM_TREE_RETENTION_SUMMARY_FILENAME
        ),
        slurm_merge_check_rows=_read_tsv_rows(
            bundle_root / _SLURM_MERGE_CHECKS_FILENAME
        ),
        slurm_merge_variant_rows=_read_tsv_rows(
            bundle_root / _SLURM_MERGE_VARIANTS_FILENAME
        ),
        slurm_merge_summary=_load_json(bundle_root / _SLURM_MERGE_SUMMARY_FILENAME),
        slurm_output_freshness_rows=_read_tsv_rows(
            bundle_root / _SLURM_OUTPUT_FRESHNESS_FILENAME
        ),
        slurm_output_freshness_check_rows=_read_tsv_rows(
            bundle_root / _SLURM_OUTPUT_FRESHNESS_CHECKS_FILENAME
        ),
        slurm_output_freshness_summary=_load_json(
            bundle_root / _SLURM_OUTPUT_FRESHNESS_SUMMARY_FILENAME
        ),
        slurm_job_status_rows=_read_tsv_rows(bundle_root / _SLURM_JOB_STATUS_FILENAME),
        slurm_partition_status_rows=_read_tsv_rows(
            bundle_root / _SLURM_PARTITION_STATUS_FILENAME
        ),
        slurm_workflow_status=_load_json(bundle_root / _SLURM_WORKFLOW_STATUS_FILENAME),
        slurm_failure_recovery_job_rows=_read_tsv_rows(
            bundle_root / _SLURM_FAILURE_RECOVERY_JOBS_FILENAME
        ),
        slurm_failure_recovery_partition_rows=_read_tsv_rows(
            bundle_root / _SLURM_FAILURE_RECOVERY_PARTITIONS_FILENAME
        ),
        slurm_failure_recovery_summary=_load_json(
            bundle_root / _SLURM_FAILURE_RECOVERY_SUMMARY_FILENAME
        ),
    )
