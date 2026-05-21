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
    load_rabies_method_sensitivity_audit_snapshot,
)
from .review_context import build_rabies_method_sensitivity_audit_review_context
from .slurm_review import record_rabies_method_sensitivity_slurm_checks

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

    review_context = build_rabies_method_sensitivity_audit_review_context(snapshot)
    config_variants = review_context.config_variants
    manifest_task_records = review_context.manifest_task_records
    parallel_summary_rows = review_context.parallel_summary_rows
    variant_summary_rows = review_context.variant_summary_rows
    logged_variants = review_context.logged_variants
    config_variant_ids = review_context.config_variant_ids
    manifest_variant_ids = review_context.manifest_variant_ids
    parallel_variant_ids = review_context.parallel_variant_ids
    summary_variant_ids = review_context.summary_variant_ids
    logged_variant_ids = review_context.logged_variant_ids
    written_variant_ids = review_context.written_variant_ids
    expected_variant_count = review_context.expected_variant_count
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
    record_rabies_method_sensitivity_slurm_checks(
        checks=checks,
        snapshot=snapshot,
        review_context=review_context,
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
