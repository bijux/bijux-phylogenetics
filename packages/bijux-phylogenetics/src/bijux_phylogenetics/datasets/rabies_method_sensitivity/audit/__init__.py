from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.common import file_sha256
from .contracts import (
    RabiesMethodSensitivityReproducibilityAuditReport,
    RabiesMethodSensitivityReproducibilityCheckRow,
    RabiesMethodSensitivityVariantAuditRow,
)
from .io import (
    _alignment_length,
    _directory_digest,
    _format_float,
    _load_json,
    _parse_task_log,
    _read_tsv_rows,
    write_rabies_method_sensitivity_reproducibility_audit_json,
    write_rabies_method_sensitivity_reproducibility_checks_table,
    write_rabies_method_sensitivity_variant_audit_table,
)
from .inventory import (
    _EXPECTED_VARIANT_FILENAMES,
    _WORKFLOW_MANIFEST_FILENAME,
    RabiesMethodSensitivityAuditSnapshot,
    load_rabies_method_sensitivity_audit_snapshot,
)

__all__ = [
    "RabiesMethodSensitivityReproducibilityAuditReport",
    "RabiesMethodSensitivityReproducibilityCheckRow",
    "RabiesMethodSensitivityVariantAuditRow",
    "audit_rabies_method_sensitivity_workflow_bundle",
    "write_rabies_method_sensitivity_reproducibility_audit_json",
    "write_rabies_method_sensitivity_reproducibility_checks_table",
    "write_rabies_method_sensitivity_variant_audit_table",
]

def audit_rabies_method_sensitivity_workflow_bundle(
    bundle_root: Path,
    *,
    sequences_path: Path,
    metadata_path: Path,
) -> RabiesMethodSensitivityReproducibilityAuditReport:
    """Audit one written bundle against current inputs, settings, and outputs."""
    snapshot = load_rabies_method_sensitivity_audit_snapshot(bundle_root)
    bundle_root = snapshot.bundle_root
    workflow_manifest_path = snapshot.workflow_manifest_path
    report_manifest_path = snapshot.report_manifest_path
    config_path = snapshot.config_path
    task_logs_root = snapshot.task_logs_root
    variants_root = snapshot.variants_root
    slurm_job_evidence_root = snapshot.slurm_job_evidence_root
    workflow_manifest = snapshot.workflow_manifest
    report_manifest = snapshot.report_manifest
    resolved_config = snapshot.resolved_config
    parallel_rows = snapshot.parallel_rows
    variant_rows = snapshot.variant_rows
    slurm_array_partition_rows = snapshot.slurm_array_partition_rows
    slurm_array_member_rows = snapshot.slurm_array_member_rows
    slurm_job_evidence_rows = snapshot.slurm_job_evidence_rows
    slurm_job_evidence_summary = snapshot.slurm_job_evidence_summary
    slurm_storage_category_rows = snapshot.slurm_storage_category_rows
    slurm_storage_variant_rows = snapshot.slurm_storage_variant_rows
    slurm_storage_summary = snapshot.slurm_storage_summary
    slurm_output_explosion_check_rows = snapshot.slurm_output_explosion_check_rows
    slurm_output_explosion_variant_rows = snapshot.slurm_output_explosion_variant_rows
    slurm_output_explosion_summary = snapshot.slurm_output_explosion_summary
    slurm_tree_retention_check_rows = snapshot.slurm_tree_retention_check_rows
    slurm_tree_retention_file_rows = snapshot.slurm_tree_retention_file_rows
    slurm_tree_retention_summary = snapshot.slurm_tree_retention_summary
    slurm_merge_check_rows = snapshot.slurm_merge_check_rows
    slurm_merge_variant_rows = snapshot.slurm_merge_variant_rows
    slurm_merge_summary = snapshot.slurm_merge_summary
    slurm_output_freshness_rows = snapshot.slurm_output_freshness_rows
    slurm_output_freshness_check_rows = snapshot.slurm_output_freshness_check_rows
    slurm_output_freshness_summary = snapshot.slurm_output_freshness_summary
    slurm_job_status_rows = snapshot.slurm_job_status_rows
    slurm_partition_status_rows = snapshot.slurm_partition_status_rows
    slurm_workflow_status = snapshot.slurm_workflow_status
    slurm_failure_recovery_job_rows = snapshot.slurm_failure_recovery_job_rows
    slurm_failure_recovery_partition_rows = snapshot.slurm_failure_recovery_partition_rows
    slurm_failure_recovery_summary = snapshot.slurm_failure_recovery_summary

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
