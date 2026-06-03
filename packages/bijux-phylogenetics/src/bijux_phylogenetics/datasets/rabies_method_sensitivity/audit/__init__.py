from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.common import file_sha256

from .contracts import (
    RabiesMethodSensitivityReproducibilityAuditReport,
    RabiesMethodSensitivityReproducibilityCheckRow,
    RabiesMethodSensitivityVariantAuditRow,
)
from .inventory import (
    _WORKFLOW_MANIFEST_FILENAME,
    load_rabies_method_sensitivity_audit_snapshot,
)
from .io import (
    _load_json,
    write_rabies_method_sensitivity_reproducibility_audit_json,
    write_rabies_method_sensitivity_reproducibility_checks_table,
    write_rabies_method_sensitivity_variant_audit_table,
)
from .review_context import build_rabies_method_sensitivity_audit_review_context
from .slurm_review import record_rabies_method_sensitivity_slurm_checks
from .variant_review import build_rabies_method_sensitivity_variant_audit_rows

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
    workflow_manifest = snapshot.workflow_manifest
    report_manifest = snapshot.report_manifest
    resolved_config = snapshot.resolved_config
    parallel_rows = snapshot.parallel_rows
    variant_rows = snapshot.variant_rows

    checks: list[RabiesMethodSensitivityReproducibilityCheckRow] = []

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

    variant_audit_rows = build_rabies_method_sensitivity_variant_audit_rows(
        snapshot=snapshot,
        review_context=review_context,
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
