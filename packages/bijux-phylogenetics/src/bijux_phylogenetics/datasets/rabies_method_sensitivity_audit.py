from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.engines.common import file_sha256
from bijux_phylogenetics.io.fasta._shared import load_fasta_alignment

__all__ = [
    "RabiesMethodSensitivityReproducibilityAuditReport",
    "RabiesMethodSensitivityReproducibilityCheckRow",
    "RabiesMethodSensitivityVariantAuditRow",
    "audit_rabies_method_sensitivity_workflow_bundle",
    "write_rabies_method_sensitivity_reproducibility_audit_json",
    "write_rabies_method_sensitivity_reproducibility_checks_table",
    "write_rabies_method_sensitivity_variant_audit_table",
]

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
class RabiesMethodSensitivityReproducibilityCheckRow:
    """One machine-readable pass/fail check within the bundle audit."""

    check_id: str
    surface: str
    status: str
    expected: str
    observed: str
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityVariantAuditRow:
    """One per-variant provenance and file-inventory summary."""

    variant_id: str
    status: str
    output_file_count: int
    output_byte_count: int
    output_digest: str
    missing_required_files: tuple[str, ...]
    unexpected_files: tuple[str, ...]
    issues: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityReproducibilityAuditReport:
    """One reviewer-facing reproducibility audit for the workflow bundle."""

    dataset_id: str
    bundle_root: Path
    workflow_manifest_path: Path
    report_manifest_path: Path
    config_path: Path
    sequences_path: Path
    metadata_path: Path
    all_passed: bool
    check_count: int
    failed_check_count: int
    variant_count: int
    failed_variant_count: int
    checks: tuple[RabiesMethodSensitivityReproducibilityCheckRow, ...]
    variants: tuple[RabiesMethodSensitivityVariantAuditRow, ...]


def audit_rabies_method_sensitivity_workflow_bundle(
    bundle_root: Path,
    *,
    sequences_path: Path,
    metadata_path: Path,
) -> RabiesMethodSensitivityReproducibilityAuditReport:
    """Audit one written bundle against current inputs, settings, and outputs."""
    bundle_root = bundle_root.resolve()
    workflow_manifest_path = bundle_root / _WORKFLOW_MANIFEST_FILENAME
    report_manifest_path = bundle_root / _REPORT_MANIFEST_FILENAME
    config_path = bundle_root / _CONFIG_FILENAME
    parallel_summary_path = bundle_root / _PARALLEL_SUMMARY_FILENAME
    variant_summary_path = bundle_root / _VARIANT_SUMMARY_FILENAME
    slurm_array_partitions_path = bundle_root / _SLURM_ARRAY_PARTITIONS_FILENAME
    slurm_array_members_path = bundle_root / _SLURM_ARRAY_MEMBERS_FILENAME
    slurm_job_evidence_root = bundle_root / _SLURM_JOB_EVIDENCE_DIRECTORY
    slurm_job_evidence_index_path = bundle_root / _SLURM_JOB_EVIDENCE_INDEX_FILENAME
    slurm_job_evidence_summary_path = (
        bundle_root / _SLURM_JOB_EVIDENCE_SUMMARY_FILENAME
    )
    slurm_storage_categories_path = bundle_root / _SLURM_STORAGE_CATEGORIES_FILENAME
    slurm_storage_variants_path = bundle_root / _SLURM_STORAGE_VARIANTS_FILENAME
    slurm_storage_summary_path = bundle_root / _SLURM_STORAGE_SUMMARY_FILENAME
    slurm_output_explosion_checks_path = (
        bundle_root / _SLURM_OUTPUT_EXPLOSION_CHECKS_FILENAME
    )
    slurm_output_explosion_variants_path = (
        bundle_root / _SLURM_OUTPUT_EXPLOSION_VARIANTS_FILENAME
    )
    slurm_output_explosion_summary_path = (
        bundle_root / _SLURM_OUTPUT_EXPLOSION_SUMMARY_FILENAME
    )
    slurm_tree_retention_checks_path = (
        bundle_root / _SLURM_TREE_RETENTION_CHECKS_FILENAME
    )
    slurm_tree_retention_files_path = bundle_root / _SLURM_TREE_RETENTION_FILES_FILENAME
    slurm_tree_retention_summary_path = (
        bundle_root / _SLURM_TREE_RETENTION_SUMMARY_FILENAME
    )
    slurm_merge_checks_path = bundle_root / _SLURM_MERGE_CHECKS_FILENAME
    slurm_merge_variants_path = bundle_root / _SLURM_MERGE_VARIANTS_FILENAME
    slurm_merge_summary_path = bundle_root / _SLURM_MERGE_SUMMARY_FILENAME
    slurm_output_freshness_path = bundle_root / _SLURM_OUTPUT_FRESHNESS_FILENAME
    slurm_output_freshness_checks_path = (
        bundle_root / _SLURM_OUTPUT_FRESHNESS_CHECKS_FILENAME
    )
    slurm_output_freshness_summary_path = (
        bundle_root / _SLURM_OUTPUT_FRESHNESS_SUMMARY_FILENAME
    )
    slurm_job_status_path = bundle_root / _SLURM_JOB_STATUS_FILENAME
    slurm_partition_status_path = bundle_root / _SLURM_PARTITION_STATUS_FILENAME
    slurm_workflow_status_path = bundle_root / _SLURM_WORKFLOW_STATUS_FILENAME
    slurm_failure_recovery_jobs_path = bundle_root / _SLURM_FAILURE_RECOVERY_JOBS_FILENAME
    slurm_failure_recovery_partitions_path = (
        bundle_root / _SLURM_FAILURE_RECOVERY_PARTITIONS_FILENAME
    )
    slurm_failure_recovery_summary_path = (
        bundle_root / _SLURM_FAILURE_RECOVERY_SUMMARY_FILENAME
    )
    task_logs_root = bundle_root / _TASK_LOGS_DIRECTORY
    variants_root = bundle_root / _VARIANTS_DIRECTORY

    workflow_manifest = _load_json(workflow_manifest_path)
    report_manifest = _load_json(report_manifest_path)
    resolved_config = _load_json(config_path)
    parallel_rows = _read_tsv_rows(parallel_summary_path)
    variant_rows = _read_tsv_rows(variant_summary_path)
    slurm_array_partition_rows = _read_tsv_rows(slurm_array_partitions_path)
    slurm_array_member_rows = _read_tsv_rows(slurm_array_members_path)
    slurm_job_evidence_rows = _read_tsv_rows(slurm_job_evidence_index_path)
    slurm_job_evidence_summary = _load_json(slurm_job_evidence_summary_path)
    slurm_storage_category_rows = _read_tsv_rows(slurm_storage_categories_path)
    slurm_storage_variant_rows = _read_tsv_rows(slurm_storage_variants_path)
    slurm_storage_summary = _load_json(slurm_storage_summary_path)
    slurm_output_explosion_check_rows = _read_tsv_rows(
        slurm_output_explosion_checks_path
    )
    slurm_output_explosion_variant_rows = _read_tsv_rows(
        slurm_output_explosion_variants_path
    )
    slurm_output_explosion_summary = _load_json(slurm_output_explosion_summary_path)
    slurm_tree_retention_check_rows = _read_tsv_rows(slurm_tree_retention_checks_path)
    slurm_tree_retention_file_rows = _read_tsv_rows(slurm_tree_retention_files_path)
    slurm_tree_retention_summary = _load_json(slurm_tree_retention_summary_path)
    slurm_merge_check_rows = _read_tsv_rows(slurm_merge_checks_path)
    slurm_merge_variant_rows = _read_tsv_rows(slurm_merge_variants_path)
    slurm_merge_summary = _load_json(slurm_merge_summary_path)
    slurm_output_freshness_rows = _read_tsv_rows(slurm_output_freshness_path)
    slurm_output_freshness_check_rows = _read_tsv_rows(
        slurm_output_freshness_checks_path
    )
    slurm_output_freshness_summary = _load_json(
        slurm_output_freshness_summary_path
    )
    slurm_job_status_rows = _read_tsv_rows(slurm_job_status_path)
    slurm_partition_status_rows = _read_tsv_rows(slurm_partition_status_path)
    slurm_workflow_status = _load_json(slurm_workflow_status_path)
    slurm_failure_recovery_job_rows = _read_tsv_rows(slurm_failure_recovery_jobs_path)
    slurm_failure_recovery_partition_rows = _read_tsv_rows(
        slurm_failure_recovery_partitions_path
    )
    slurm_failure_recovery_summary = _load_json(slurm_failure_recovery_summary_path)

    checks: list[RabiesMethodSensitivityReproducibilityCheckRow] = []
    variant_audit_rows: list[RabiesMethodSensitivityVariantAuditRow] = []

    def add_check(
        check_id: str,
        *,
        surface: str,
        condition: bool,
        expected: object,
        observed: object,
        detail: str,
    ) -> None:
        checks.append(
            RabiesMethodSensitivityReproducibilityCheckRow(
                check_id=check_id,
                surface=surface,
                status="passed" if condition else "failed",
                expected="" if expected is None else str(expected),
                observed="" if observed is None else str(observed),
                detail=detail,
            )
        )

    input_paths = {
        "sequences.fasta": sequences_path,
        "metadata.csv": metadata_path,
    }
    recorded_input_checksums = {
        str(key): str(value)
        for key, value in dict(resolved_config.get("input_checksums", {})).items()
    }
    for filename, input_path in input_paths.items():
        recorded_checksum = recorded_input_checksums.get(filename)
        observed_checksum = file_sha256(input_path)
        add_check(
            f"input-checksum:{filename}",
            surface="input-checksum",
            condition=recorded_checksum == observed_checksum,
            expected=recorded_checksum,
            observed=observed_checksum,
            detail=f"{filename} matches the checksum recorded in workflow-config.resolved.json",
        )

    manifest_output_paths = {
        str(key): workflow_manifest_path.parent / Path(value)
        for key, value in dict(workflow_manifest.get("output_paths", {})).items()
    }
    manifest_output_checksums = {
        str(key): str(value)
        for key, value in dict(workflow_manifest.get("output_checksums", {})).items()
    }
    for key, output_path in sorted(manifest_output_paths.items()):
        if key in {
            "task_logs_root",
            "variants_root",
            "slurm_array_scripts_root",
            "slurm_job_evidence_root",
        }:
            add_check(
                f"workflow-manifest:{key}",
                surface="workflow-manifest",
                condition=output_path.is_dir(),
                expected="directory exists",
                observed="present" if output_path.is_dir() else "missing",
                detail=f"{key} resolves to an existing directory",
            )
            continue
        recorded_checksum = manifest_output_checksums.get(key)
        output_exists = output_path.is_file()
        observed_checksum = None if not output_exists else file_sha256(output_path)
        add_check(
            f"workflow-manifest:{key}",
            surface="workflow-manifest",
            condition=output_exists and recorded_checksum == observed_checksum,
            expected=recorded_checksum,
            observed=observed_checksum,
            detail=f"{key} matches the checksum recorded in {_WORKFLOW_MANIFEST_FILENAME}",
        )

    linked_artifacts = dict(report_manifest.get("linked_artifacts", {}))
    linked_artifact_count = int(report_manifest.get("linked_artifact_count", 0))
    add_check(
        "report-manifest:linked-artifact-count",
        surface="report-manifest",
        condition=linked_artifact_count == len(linked_artifacts),
        expected=linked_artifact_count,
        observed=len(linked_artifacts),
        detail="linked artifact count matches the linked_artifacts payload",
    )
    for key, payload in sorted(linked_artifacts.items()):
        artifact_path = (report_manifest_path.parent / str(payload["path"])).resolve()
        output_exists = artifact_path.is_file()
        observed_checksum = None if not output_exists else file_sha256(artifact_path)
        add_check(
            f"report-manifest:{key}",
            surface="report-manifest",
            condition=output_exists and str(payload["sha256"]) == observed_checksum,
            expected=payload["sha256"],
            observed=observed_checksum,
            detail=f"{key} matches the checksum recorded in the report manifest",
        )

    config_variants = {
        str(row["variant_id"]): row for row in list(resolved_config.get("variants", []))
    }
    manifest_task_records = {
        str(row["variant_id"]): row for row in list(workflow_manifest.get("task_records", []))
    }
    parallel_summary_rows = {
        str(row["variant_id"]): row for row in parallel_rows
    }
    variant_summary_rows = {
        str(row["variant_id"]): row for row in variant_rows
    }
    logged_variants = {
        path.stem: _parse_task_log(path) for path in sorted(task_logs_root.glob("*.log"))
    }

    config_variant_ids = sorted(config_variants)
    manifest_variant_ids = sorted(manifest_task_records)
    parallel_variant_ids = sorted(parallel_summary_rows)
    summary_variant_ids = sorted(variant_summary_rows)
    logged_variant_ids = sorted(logged_variants)
    written_variant_ids = sorted(
        path.name for path in variants_root.iterdir() if path.is_dir()
    )
    expected_variant_count = len(config_variant_ids)
    add_check(
        "variant-sets:manifest",
        surface="variant-sets",
        condition=config_variant_ids == manifest_variant_ids,
        expected=config_variant_ids,
        observed=manifest_variant_ids,
        detail="workflow manifest task records cover the configured variant ids",
    )
    add_check(
        "variant-sets:parallel-summary",
        surface="variant-sets",
        condition=config_variant_ids == parallel_variant_ids,
        expected=config_variant_ids,
        observed=parallel_variant_ids,
        detail="parallel summary rows cover the configured variant ids",
    )
    add_check(
        "variant-sets:variant-summary",
        surface="variant-sets",
        condition=config_variant_ids == summary_variant_ids,
        expected=config_variant_ids,
        observed=summary_variant_ids,
        detail="variant summary rows cover the configured variant ids",
    )
    add_check(
        "variant-sets:task-logs",
        surface="variant-sets",
        condition=config_variant_ids == logged_variant_ids,
        expected=config_variant_ids,
        observed=logged_variant_ids,
        detail="task log files cover the configured variant ids",
    )
    add_check(
        "variant-sets:variant-directories",
        surface="variant-sets",
        condition=config_variant_ids == written_variant_ids,
        expected=config_variant_ids,
        observed=written_variant_ids,
        detail="variant output directories cover the configured variant ids",
    )
    add_check(
        "variant-count:parallel-summary",
        surface="variant-count",
        condition=expected_variant_count == len(parallel_rows),
        expected=expected_variant_count,
        observed=len(parallel_rows),
        detail="parallel summary row count matches the configured variant count",
    )
    add_check(
        "variant-count:variant-summary",
        surface="variant-count",
        condition=expected_variant_count == len(variant_rows),
        expected=expected_variant_count,
        observed=len(variant_rows),
        detail="variant summary row count matches the configured variant count",
    )
    slurm_partition_ids = {
        str(row["partition_id"]) for row in slurm_array_partition_rows
    }
    slurm_script_paths = {
        str(row["script_path"]) for row in slurm_array_partition_rows
    }
    member_partition_ids = {
        str(row["partition_id"]) for row in slurm_array_member_rows
    }
    member_variant_ids = sorted(str(row["variant_id"]) for row in slurm_array_member_rows)
    add_check(
        "slurm-arrays:partition-coverage",
        surface="slurm-arrays",
        condition=bool(slurm_partition_ids) and slurm_partition_ids == member_partition_ids,
        expected=sorted(slurm_partition_ids),
        observed=sorted(member_partition_ids),
        detail="array member rows cover the same partition ids as the partition summary",
    )
    add_check(
        "slurm-arrays:member-coverage",
        surface="slurm-arrays",
        condition=config_variant_ids == member_variant_ids,
        expected=config_variant_ids,
        observed=member_variant_ids,
        detail="array member rows cover the configured variant ids",
    )
    slurm_job_evidence_variant_ids = sorted(
        str(row["variant_id"]) for row in slurm_job_evidence_rows
    )
    add_check(
        "slurm-job-evidence:job-coverage",
        surface="slurm-job-evidence",
        condition=config_variant_ids == slurm_job_evidence_variant_ids,
        expected=config_variant_ids,
        observed=slurm_job_evidence_variant_ids,
        detail="job-evidence rows cover the configured variant ids",
    )
    add_check(
        "slurm-job-evidence:job-count",
        surface="slurm-job-evidence",
        condition=int(slurm_job_evidence_summary["job_count"])
        == len(slurm_job_evidence_rows),
        expected=slurm_job_evidence_summary["job_count"],
        observed=len(slurm_job_evidence_rows),
        detail="job-evidence summary job_count matches the written job-evidence rows",
    )
    job_evidence_artifact_file_count = len(
        [path for path in slurm_job_evidence_root.rglob("*") if path.is_file()]
    )
    add_check(
        "slurm-job-evidence:artifact-file-count",
        surface="slurm-job-evidence",
        condition=int(slurm_job_evidence_summary["total_artifact_file_count"])
        == job_evidence_artifact_file_count,
        expected=slurm_job_evidence_summary["total_artifact_file_count"],
        observed=job_evidence_artifact_file_count,
        detail="job-evidence summary total_artifact_file_count matches the written evidence files",
    )
    slurm_storage_category_ids = sorted(
        str(row["category_id"]) for row in slurm_storage_category_rows
    )
    slurm_storage_variant_ids = sorted(
        str(row["variant_id"]) for row in slurm_storage_variant_rows
    )
    add_check(
        "slurm-storage:category-coverage",
        surface="slurm-storage",
        condition=slurm_storage_category_ids
        == ["logs", "outputs", "posterior_samples", "reports", "trees"],
        expected=["logs", "outputs", "posterior_samples", "reports", "trees"],
        observed=slurm_storage_category_ids,
        detail="storage category rows cover the explicit retained-storage categories",
    )
    add_check(
        "slurm-storage:variant-coverage",
        surface="slurm-storage",
        condition=config_variant_ids == slurm_storage_variant_ids,
        expected=config_variant_ids,
        observed=slurm_storage_variant_ids,
        detail="storage variant rows cover the configured variant ids",
    )
    storage_total_file_count = sum(
        int(row["total_file_count"]) for row in slurm_storage_category_rows
    )
    storage_total_byte_count = sum(
        int(row["total_byte_count"]) for row in slurm_storage_category_rows
    )
    add_check(
        "slurm-storage:summary-totals",
        surface="slurm-storage",
        condition=int(slurm_storage_summary["total_file_count"])
        == storage_total_file_count
        and int(slurm_storage_summary["total_byte_count"]) == storage_total_byte_count,
        expected=(
            int(slurm_storage_summary["total_file_count"]),
            int(slurm_storage_summary["total_byte_count"]),
        ),
        observed=(storage_total_file_count, storage_total_byte_count),
        detail="storage summary totals match the written category ledger",
    )
    largest_storage_variant_row = max(
        slurm_storage_variant_rows,
        key=lambda row: (int(row["total_byte_count"]), str(row["variant_id"])),
    )
    add_check(
        "slurm-storage:largest-variant",
        surface="slurm-storage",
        condition=str(slurm_storage_summary["largest_variant_id"])
        == str(largest_storage_variant_row["variant_id"])
        and int(slurm_storage_summary["largest_variant_total_byte_count"])
        == int(largest_storage_variant_row["total_byte_count"]),
        expected=(
            str(slurm_storage_summary["largest_variant_id"]),
            int(slurm_storage_summary["largest_variant_total_byte_count"]),
        ),
        observed=(
            str(largest_storage_variant_row["variant_id"]),
            int(largest_storage_variant_row["total_byte_count"]),
        ),
        detail="storage summary largest-variant fields match the written per-variant storage ledger",
    )
    slurm_output_explosion_variant_ids = sorted(
        str(row["variant_id"]) for row in slurm_output_explosion_variant_rows
    )
    add_check(
        "slurm-output-explosion:variant-coverage",
        surface="slurm-output-explosion",
        condition=config_variant_ids == slurm_output_explosion_variant_ids,
        expected=config_variant_ids,
        observed=slurm_output_explosion_variant_ids,
        detail="output-explosion rows cover the configured variant ids",
    )
    add_check(
        "slurm-output-explosion:check-count",
        surface="slurm-output-explosion",
        condition=int(slurm_output_explosion_summary["check_count"])
        == len(slurm_output_explosion_check_rows),
        expected=slurm_output_explosion_summary["check_count"],
        observed=len(slurm_output_explosion_check_rows),
        detail="output-explosion summary check_count matches the written check rows",
    )
    observed_output_explosion_variant_total = (
        sum(
            1
            for row in slurm_output_explosion_variant_rows
            if str(row["risk_status"]) == "low"
        )
        + sum(
            1
            for row in slurm_output_explosion_variant_rows
            if str(row["risk_status"]) == "warning"
        )
        + sum(
            1
            for row in slurm_output_explosion_variant_rows
            if str(row["risk_status"]) == "high"
        )
    )
    add_check(
        "slurm-output-explosion:variant-counts",
        surface="slurm-output-explosion",
        condition=int(slurm_output_explosion_summary["variant_count"])
        == observed_output_explosion_variant_total,
        expected=slurm_output_explosion_summary["variant_count"],
        observed=observed_output_explosion_variant_total,
        detail="output-explosion summary variant_count matches the risk-classified variant rows",
    )
    observed_overall_output_explosion_status = "high"
    if (
        int(slurm_output_explosion_summary["failed_check_count"]) == 0
        and int(slurm_output_explosion_summary["high_risk_variant_count"]) == 0
        and int(slurm_output_explosion_summary["global_issue_count"]) == 0
        and int(slurm_output_explosion_summary["warning_variant_count"]) == 0
    ):
        observed_overall_output_explosion_status = "low"
    elif (
        int(slurm_output_explosion_summary["failed_check_count"]) == 0
        and int(slurm_output_explosion_summary["high_risk_variant_count"]) == 0
    ):
        observed_overall_output_explosion_status = "warning"
    add_check(
        "slurm-output-explosion:overall-risk",
        surface="slurm-output-explosion",
        condition=str(slurm_output_explosion_summary["overall_risk_status"])
        == observed_overall_output_explosion_status,
        expected=slurm_output_explosion_summary["overall_risk_status"],
        observed=observed_overall_output_explosion_status,
        detail="output-explosion summary overall_risk_status matches the written risk counts",
    )
    slurm_tree_retention_variant_ids = sorted(
        {str(row["variant_id"]) for row in slurm_tree_retention_file_rows}
    )
    add_check(
        "slurm-tree-retention:variant-coverage",
        surface="slurm-tree-retention",
        condition=config_variant_ids == slurm_tree_retention_variant_ids,
        expected=config_variant_ids,
        observed=slurm_tree_retention_variant_ids,
        detail="tree-retention file rows cover the configured variant ids",
    )
    add_check(
        "slurm-tree-retention:check-count",
        surface="slurm-tree-retention",
        condition=int(slurm_tree_retention_summary["check_count"])
        == len(slurm_tree_retention_check_rows),
        expected=slurm_tree_retention_summary["check_count"],
        observed=len(slurm_tree_retention_check_rows),
        detail="tree-retention summary check_count matches the written check rows",
    )
    observed_tree_set_file_count = sum(
        1 for row in slurm_tree_retention_file_rows if int(row["tree_count"]) > 1
    )
    add_check(
        "slurm-tree-retention:tree-set-count",
        surface="slurm-tree-retention",
        condition=int(slurm_tree_retention_summary["tree_set_file_count"])
        == observed_tree_set_file_count,
        expected=slurm_tree_retention_summary["tree_set_file_count"],
        observed=observed_tree_set_file_count,
        detail="tree-retention summary tree_set_file_count matches the inspected file rows",
    )
    observed_tree_retention_status = "required"
    if (
        int(slurm_tree_retention_summary["failed_check_count"]) == 0
        and int(slurm_tree_retention_summary["thinning_required_file_count"]) == 0
        and int(slurm_tree_retention_summary["compression_required_file_count"]) == 0
        and int(slurm_tree_retention_summary["thinning_recommended_file_count"]) == 0
        and int(slurm_tree_retention_summary["compression_recommended_file_count"]) == 0
    ):
        observed_tree_retention_status = "no_action"
    elif int(slurm_tree_retention_summary["failed_check_count"]) == 0 and (
        int(slurm_tree_retention_summary["thinning_required_file_count"]) == 0
        and int(slurm_tree_retention_summary["compression_required_file_count"]) == 0
    ):
        observed_tree_retention_status = "recommended"
    add_check(
        "slurm-tree-retention:overall-status",
        surface="slurm-tree-retention",
        condition=str(slurm_tree_retention_summary["overall_policy_status"])
        == observed_tree_retention_status,
        expected=slurm_tree_retention_summary["overall_policy_status"],
        observed=observed_tree_retention_status,
        detail="tree-retention summary overall_policy_status matches the written policy counts",
    )
    slurm_merge_variant_ids = sorted(
        str(row["variant_id"]) for row in slurm_merge_variant_rows
    )
    add_check(
        "slurm-merge:job-coverage",
        surface="slurm-merge",
        condition=config_variant_ids == slurm_merge_variant_ids,
        expected=config_variant_ids,
        observed=slurm_merge_variant_ids,
        detail="merge-variant rows cover the configured variant ids",
    )
    add_check(
        "slurm-merge:check-count",
        surface="slurm-merge",
        condition=int(slurm_merge_summary["check_count"])
        == len(slurm_merge_check_rows),
        expected=slurm_merge_summary["check_count"],
        observed=len(slurm_merge_check_rows),
        detail="merge summary check_count matches the written merge checks",
    )
    merged_variant_count = sum(
        1
        for row in slurm_merge_variant_rows
        if str(row["included_in_merge"]) == "true"
    )
    add_check(
        "slurm-merge:merged-variant-count",
        surface="slurm-merge",
        condition=int(slurm_merge_summary["merged_variant_count"])
        == merged_variant_count,
        expected=slurm_merge_summary["merged_variant_count"],
        observed=merged_variant_count,
        detail="merge summary merged_variant_count matches the written merge-variant rows",
    )
    add_check(
        "slurm-merge:merge-ready",
        surface="slurm-merge",
        condition=bool(slurm_merge_summary["merge_ready"]) == (
            int(slurm_merge_summary["failed_check_count"]) == 0
            and merged_variant_count == len(config_variant_ids)
        ),
        expected=int(slurm_merge_summary["failed_check_count"]) == 0
        and merged_variant_count == len(config_variant_ids),
        observed=slurm_merge_summary["merge_ready"],
        detail="merge summary merge_ready matches the merge-check and merged-variant totals",
    )
    slurm_output_freshness_variant_ids = sorted(
        str(row["variant_id"]) for row in slurm_output_freshness_rows
    )
    add_check(
        "slurm-freshness:job-coverage",
        surface="slurm-freshness",
        condition=config_variant_ids == slurm_output_freshness_variant_ids,
        expected=config_variant_ids,
        observed=slurm_output_freshness_variant_ids,
        detail="output-freshness rows cover the configured variant ids",
    )
    freshness_failed_check_count = sum(
        1
        for row in slurm_output_freshness_check_rows
        if str(row["status"]) == "failed"
    )
    freshness_stale_job_count = sum(
        1
        for row in slurm_output_freshness_rows
        if str(row["freshness_status"]) == "stale"
    )
    add_check(
        "slurm-freshness:check-count",
        surface="slurm-freshness",
        condition=int(slurm_output_freshness_summary["check_count"])
        == len(slurm_output_freshness_check_rows),
        expected=slurm_output_freshness_summary["check_count"],
        observed=len(slurm_output_freshness_check_rows),
        detail="output-freshness summary check_count matches the written check rows",
    )
    add_check(
        "slurm-freshness:job-count",
        surface="slurm-freshness",
        condition=(
            int(slurm_output_freshness_summary["fresh_job_count"])
            + int(slurm_output_freshness_summary["stale_job_count"])
        )
        == len(slurm_output_freshness_rows),
        expected=len(slurm_output_freshness_rows),
        observed=(
            int(slurm_output_freshness_summary["fresh_job_count"])
            + int(slurm_output_freshness_summary["stale_job_count"])
        ),
        detail="output-freshness summary job counts match the written per-job ledger",
    )
    add_check(
        "slurm-freshness:all-outputs-fresh",
        surface="slurm-freshness",
        condition=bool(slurm_output_freshness_summary["all_outputs_fresh"]) == (
            freshness_failed_check_count == 0 and freshness_stale_job_count == 0
        ),
        expected=freshness_failed_check_count == 0 and freshness_stale_job_count == 0,
        observed=slurm_output_freshness_summary["all_outputs_fresh"],
        detail="output-freshness summary all_outputs_fresh matches the written checks and per-job statuses",
    )
    for script_path_text in sorted(slurm_script_paths):
        script_path = bundle_root / script_path_text
        add_check(
            f"slurm-arrays:script:{Path(script_path_text).name}",
            surface="slurm-arrays",
            condition=script_path.is_file(),
            expected="script file exists",
            observed="present" if script_path.is_file() else "missing",
            detail="array partition table references an existing sbatch script",
        )
    slurm_job_status_variant_ids = sorted(
        str(row["variant_id"]) for row in slurm_job_status_rows
    )
    slurm_partition_status_ids = sorted(
        str(row["partition_id"]) for row in slurm_partition_status_rows
    )
    slurm_output_freshness_rows_by_variant = {
        str(row["variant_id"]): row for row in slurm_output_freshness_rows
    }
    slurm_job_evidence_rows_by_variant = {
        str(row["variant_id"]): row for row in slurm_job_evidence_rows
    }
    slurm_merge_rows_by_variant = {
        str(row["variant_id"]): row for row in slurm_merge_variant_rows
    }
    slurm_job_status_rows_by_variant = {
        str(row["variant_id"]): row for row in slurm_job_status_rows
    }
    add_check(
        "slurm-status:job-coverage",
        surface="slurm-status",
        condition=config_variant_ids == slurm_job_status_variant_ids,
        expected=config_variant_ids,
        observed=slurm_job_status_variant_ids,
        detail="job-status rows cover the configured variant ids",
    )
    add_check(
        "slurm-status:partition-coverage",
        surface="slurm-status",
        condition=sorted(slurm_partition_ids) == slurm_partition_status_ids,
        expected=sorted(slurm_partition_ids),
        observed=slurm_partition_status_ids,
        detail="partition-status rows cover the same partition ids as the array partition table",
    )
    add_check(
        "slurm-status:workflow-job-count",
        surface="slurm-status",
        condition=int(slurm_workflow_status["job_count"]) == len(slurm_job_status_rows),
        expected=slurm_workflow_status["job_count"],
        observed=len(slurm_job_status_rows),
        detail="workflow status job_count matches the number of job-status rows",
    )
    add_check(
        "slurm-status:workflow-partition-count",
        surface="slurm-status",
        condition=int(slurm_workflow_status["partition_count"])
        == len(slurm_partition_status_rows),
        expected=slurm_workflow_status["partition_count"],
        observed=len(slurm_partition_status_rows),
        detail="workflow status partition_count matches the number of partition-status rows",
    )
    slurm_failure_recovery_variant_ids = sorted(
        {str(row["variant_id"]) for row in slurm_failure_recovery_job_rows}
    )
    add_check(
        "slurm-failure-recovery:variant-coverage",
        surface="slurm-failure-recovery",
        condition=config_variant_ids == slurm_failure_recovery_variant_ids,
        expected=config_variant_ids,
        observed=slurm_failure_recovery_variant_ids,
        detail="failure-recovery job rows cover the configured variant ids",
    )
    add_check(
        "slurm-failure-recovery:partition-coverage",
        surface="slurm-failure-recovery",
        condition=int(slurm_failure_recovery_summary["partition_count"])
        == len(slurm_failure_recovery_partition_rows),
        expected=slurm_failure_recovery_summary["partition_count"],
        observed=len(slurm_failure_recovery_partition_rows),
        detail="failure-recovery summary partition_count matches the written partition rows",
    )
    observed_rerunnable_job_count = sum(
        1 for row in slurm_failure_recovery_job_rows if str(row["rerunnable"]) == "true"
    )
    add_check(
        "slurm-failure-recovery:rerunnable-count",
        surface="slurm-failure-recovery",
        condition=int(slurm_failure_recovery_summary["rerunnable_job_count"])
        == observed_rerunnable_job_count,
        expected=slurm_failure_recovery_summary["rerunnable_job_count"],
        observed=observed_rerunnable_job_count,
        detail="failure-recovery summary rerunnable_job_count matches the written job rows",
    )
    observed_failure_recovery_status = "clean"
    if observed_rerunnable_job_count > 0:
        observed_failure_recovery_status = "recovery_needed"
    elif int(slurm_failure_recovery_summary["blocked_job_count"]) > 0:
        observed_failure_recovery_status = "workflow_active"
    add_check(
        "slurm-failure-recovery:overall-status",
        surface="slurm-failure-recovery",
        condition=str(slurm_failure_recovery_summary["overall_recovery_status"])
        == observed_failure_recovery_status,
        expected=slurm_failure_recovery_summary["overall_recovery_status"],
        observed=observed_failure_recovery_status,
        detail="failure-recovery summary overall_recovery_status matches the written rerun decisions",
    )
    for variant_id in config_variant_ids:
        freshness_row = slurm_output_freshness_rows_by_variant.get(variant_id)
        job_evidence_row = slurm_job_evidence_rows_by_variant.get(variant_id)
        merge_row = slurm_merge_rows_by_variant.get(variant_id)
        job_status_row = slurm_job_status_rows_by_variant.get(variant_id)
        add_check(
            f"slurm-freshness:status-link:{variant_id}",
            surface="slurm-freshness",
            condition=(
                freshness_row is not None
                and job_status_row is not None
                and str(freshness_row["freshness_status"])
                == str(job_status_row["output_freshness_status"])
            ),
            expected=(
                None
                if freshness_row is None
                else freshness_row["freshness_status"]
            ),
            observed=(
                None
                if job_status_row is None
                else job_status_row["output_freshness_status"]
            ),
            detail="job-status rows expose the same output-freshness status as the freshness ledger",
        )
        add_check(
            f"slurm-job-evidence:status-link:{variant_id}",
            surface="slurm-job-evidence",
            condition=(
                job_evidence_row is not None
                and job_status_row is not None
                and str(job_evidence_row["status"])
                == str(job_status_row["task_status"])
            ),
            expected=(
                None if job_status_row is None else job_status_row["task_status"]
            ),
            observed=(
                None if job_evidence_row is None else job_evidence_row["status"]
            ),
            detail="job-evidence rows expose the same terminal task status recorded in the job-status ledger",
        )
        add_check(
            f"slurm-merge:status-link:{variant_id}",
            surface="slurm-merge",
            condition=(
                merge_row is not None
                and job_status_row is not None
                and str(merge_row["job_status"]) == str(job_status_row["status"])
            ),
            expected=None if job_status_row is None else job_status_row["status"],
            observed=None if merge_row is None else merge_row["job_status"],
            detail="merge-variant rows expose the same job status as the job-status ledger",
        )
        add_check(
            f"slurm-merge:freshness-link:{variant_id}",
            surface="slurm-merge",
            condition=(
                merge_row is not None
                and freshness_row is not None
                and str(merge_row["output_freshness_status"])
                == str(freshness_row["freshness_status"])
            ),
            expected=(
                None if freshness_row is None else freshness_row["freshness_status"]
            ),
            observed=(
                None if merge_row is None else merge_row["output_freshness_status"]
            ),
            detail="merge-variant rows expose the same freshness status as the freshness ledger",
        )
        if job_evidence_row is None:
            continue
        evidence_json_path = bundle_root / str(job_evidence_row["evidence_json_path"])
        evidence_html_path = bundle_root / str(job_evidence_row["evidence_html_path"])
        add_check(
            f"slurm-job-evidence:json:{variant_id}",
            surface="slurm-job-evidence",
            condition=evidence_json_path.is_file(),
            expected="file exists",
            observed="present" if evidence_json_path.is_file() else "missing",
            detail="job-evidence index references an existing JSON evidence package",
        )
        add_check(
            f"slurm-job-evidence:html:{variant_id}",
            surface="slurm-job-evidence",
            condition=evidence_html_path.is_file(),
            expected="file exists",
            observed="present" if evidence_html_path.is_file() else "missing",
            detail="job-evidence index references an existing HTML evidence package",
        )

    for variant_id in config_variant_ids:
        config_row = config_variants[variant_id]
        manifest_row = manifest_task_records.get(variant_id)
        parallel_row = parallel_summary_rows.get(variant_id)
        summary_row = variant_summary_rows.get(variant_id)
        log_row = logged_variants.get(variant_id)
        variant_root = variants_root / variant_id
        issues: list[str] = []

        if manifest_row is None:
            issues.append("workflow manifest task record is missing")
        if parallel_row is None:
            issues.append("parallel execution summary row is missing")
        if summary_row is None:
            issues.append("variant summary row is missing")
        if log_row is None:
            issues.append("task log is missing")

        expected_output_root = Path("variants", variant_id).as_posix()
        if manifest_row is not None and str(manifest_row.get("output_root")) != expected_output_root:
            issues.append("workflow manifest output_root differs from the expected variant directory")
        if parallel_row is not None and str(parallel_row.get("log_path")) != Path("parallel-logs", f"{variant_id}.log").as_posix():
            issues.append("parallel summary log_path differs from the expected task log path")
        if log_row is not None and str(log_row.get("output_root")) != expected_output_root:
            issues.append("task log output_root differs from the expected variant directory")
        if manifest_row is not None and str(manifest_row.get("status")) != "succeeded":
            issues.append("workflow manifest does not record the variant as succeeded")
        if parallel_row is not None and str(parallel_row.get("status")) != "succeeded":
            issues.append("parallel execution summary does not record the variant as succeeded")
        if log_row is not None and str(log_row.get("status")) != "succeeded":
            issues.append("task log does not record the variant as succeeded")

        for field_name in ("alignment_mode", "trimming_mode"):
            expected_value = str(config_row[field_name])
            if log_row is not None and str(log_row.get(field_name)) != expected_value:
                issues.append(f"task log {field_name} does not match the resolved config")
            if summary_row is not None and str(summary_row.get(field_name)) != expected_value:
                issues.append(f"variant summary {field_name} does not match the resolved config")

        expected_trim_gap_threshold = _format_float(float(config_row["trim_gap_threshold"]))
        if log_row is not None and str(log_row.get("trim_gap_threshold")) != expected_trim_gap_threshold:
            issues.append("task log trim_gap_threshold does not match the resolved config")
        if summary_row is not None and str(summary_row.get("trim_gap_threshold")) != expected_trim_gap_threshold:
            issues.append(
                "variant summary trim_gap_threshold does not match the resolved config"
            )

        variant_output_paths = sorted(
            path for path in variant_root.iterdir() if path.is_file()
        ) if variant_root.is_dir() else []
        variant_filenames = tuple(path.name for path in variant_output_paths)
        missing_required_files = tuple(
            name
            for name in (
                f"{variant_id}.aln",
                f"{variant_id}.trimmed.aln",
                *_EXPECTED_VARIANT_FILENAMES,
            )
            if not (variant_root / name).is_file()
        )
        unexpected_files = tuple(
            name
            for name in variant_filenames
            if name
            not in {
                f"{variant_id}.aln",
                f"{variant_id}.trimmed.aln",
                *_EXPECTED_VARIANT_FILENAMES,
            }
        )
        if not variant_root.is_dir():
            issues.append("variant output directory is missing")
        if missing_required_files:
            issues.append("one or more required variant output files are missing")
        if unexpected_files:
            issues.append("variant output directory contains unexpected files")

        alignment_path = variant_root / f"{variant_id}.aln"
        trimmed_alignment_path = variant_root / f"{variant_id}.trimmed.aln"
        if summary_row is not None and alignment_path.is_file():
            alignment_length = _alignment_length(alignment_path)
            if int(summary_row["alignment_length"]) != alignment_length:
                issues.append(
                    "variant summary alignment_length does not match the written alignment"
                )
        if summary_row is not None and trimmed_alignment_path.is_file():
            trimmed_alignment_length = _alignment_length(trimmed_alignment_path)
            if int(summary_row["trimmed_alignment_length"]) != trimmed_alignment_length:
                issues.append(
                    "variant summary trimmed_alignment_length does not match the written trimmed alignment"
                )

        output_byte_count = sum(path.stat().st_size for path in variant_output_paths)
        output_digest = _directory_digest(variant_root)
        variant_audit_rows.append(
            RabiesMethodSensitivityVariantAuditRow(
                variant_id=variant_id,
                status="passed" if not issues else "failed",
                output_file_count=len(variant_output_paths),
                output_byte_count=output_byte_count,
                output_digest=output_digest,
                missing_required_files=missing_required_files,
                unexpected_files=unexpected_files,
                issues=tuple(issues),
            )
        )

    failed_check_count = sum(1 for row in checks if row.status == "failed")
    failed_variant_count = sum(
        1 for row in variant_audit_rows if row.status == "failed"
    )
    return RabiesMethodSensitivityReproducibilityAuditReport(
        dataset_id=str(workflow_manifest["dataset_id"]),
        bundle_root=bundle_root,
        workflow_manifest_path=workflow_manifest_path,
        report_manifest_path=report_manifest_path,
        config_path=config_path,
        sequences_path=sequences_path,
        metadata_path=metadata_path,
        all_passed=failed_check_count == 0 and failed_variant_count == 0,
        check_count=len(checks),
        failed_check_count=failed_check_count,
        variant_count=len(variant_audit_rows),
        failed_variant_count=failed_variant_count,
        checks=tuple(checks),
        variants=tuple(variant_audit_rows),
    )


def write_rabies_method_sensitivity_reproducibility_checks_table(
    path: Path, report: RabiesMethodSensitivityReproducibilityAuditReport
) -> Path:
    """Write one tabular ledger of top-level audit checks."""
    return _write_tsv(
        path,
        fieldnames=("check_id", "surface", "status", "expected", "observed", "detail"),
        rows=[asdict(row) for row in report.checks],
    )


def write_rabies_method_sensitivity_variant_audit_table(
    path: Path, report: RabiesMethodSensitivityReproducibilityAuditReport
) -> Path:
    """Write one per-variant audit ledger."""
    return _write_tsv(
        path,
        fieldnames=(
            "variant_id",
            "status",
            "output_file_count",
            "output_byte_count",
            "output_digest",
            "missing_required_files",
            "unexpected_files",
            "issues",
        ),
        rows=[
            {
                **asdict(row),
                "missing_required_files": "; ".join(row.missing_required_files),
                "unexpected_files": "; ".join(row.unexpected_files),
                "issues": "; ".join(row.issues),
            }
            for row in report.variants
        ],
    )


def write_rabies_method_sensitivity_reproducibility_audit_json(
    path: Path, report: RabiesMethodSensitivityReproducibilityAuditReport
) -> Path:
    """Write one machine-readable JSON summary for the bundle audit."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    return path


def _alignment_length(path: Path) -> int:
    records = load_fasta_alignment(path)
    return len(records[0].sequence)


def _directory_digest(path: Path) -> str:
    if not path.is_dir():
        return ""
    lines = [
        f"{entry.relative_to(path).as_posix()}\t{file_sha256(entry)}"
        for entry in sorted(path.rglob("*"))
        if entry.is_file()
    ]
    payload = "\n".join(lines).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _format_float(value: float) -> str:
    return format(value, ".12g")


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _parse_task_log(path: Path) -> dict[str, str]:
    payload: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        payload[key.strip()] = value.strip()
    return payload


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(
    path: Path,
    *,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, object]],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: ""
                    if value is None
                    else value
                    for key, value in row.items()
                }
            )
    return path
