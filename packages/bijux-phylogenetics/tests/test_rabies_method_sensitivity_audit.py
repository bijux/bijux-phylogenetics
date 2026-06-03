from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.datasets import (
    audit_rabies_method_sensitivity_workflow_bundle,
    load_rabies_method_sensitivity_panel_dataset,
    write_rabies_method_sensitivity_reproducibility_audit_json,
    write_rabies_method_sensitivity_reproducibility_checks_table,
    write_rabies_method_sensitivity_variant_audit_table,
)


def test_audit_rabies_method_sensitivity_workflow_bundle_passes_on_packaged_outputs() -> (
    None
):
    dataset = load_rabies_method_sensitivity_panel_dataset()

    report = audit_rabies_method_sensitivity_workflow_bundle(
        dataset.reference_output_root,
        sequences_path=dataset.sequences_path,
        metadata_path=dataset.metadata_path,
    )

    assert report.all_passed is True
    assert report.failed_check_count == 0
    assert report.failed_variant_count == 0
    assert report.variant_count == 4
    assert report.check_count >= 66
    assert all(row.status == "passed" for row in report.variants)
    workflow_status_check = next(
        row
        for row in report.checks
        if row.check_id == "slurm-status:workflow-job-count"
    )
    assert workflow_status_check.status == "passed"
    freshness_check = next(
        row for row in report.checks if row.check_id == "slurm-freshness:job-coverage"
    )
    assert freshness_check.status == "passed"
    job_evidence_check = next(
        row
        for row in report.checks
        if row.check_id == "slurm-job-evidence:job-coverage"
    )
    assert job_evidence_check.status == "passed"
    storage_check = next(
        row
        for row in report.checks
        if row.check_id == "slurm-storage:category-coverage"
    )
    assert storage_check.status == "passed"
    output_explosion_check = next(
        row
        for row in report.checks
        if row.check_id == "slurm-output-explosion:variant-coverage"
    )
    assert output_explosion_check.status == "passed"
    tree_retention_check = next(
        row
        for row in report.checks
        if row.check_id == "slurm-tree-retention:variant-coverage"
    )
    assert tree_retention_check.status == "passed"
    merge_check = next(
        row for row in report.checks if row.check_id == "slurm-merge:job-coverage"
    )
    assert merge_check.status == "passed"
    failure_recovery_check = next(
        row
        for row in report.checks
        if row.check_id == "slurm-failure-recovery:variant-coverage"
    )
    assert failure_recovery_check.status == "passed"


def test_audit_rabies_method_sensitivity_workflow_bundle_detects_input_checksum_drift(
    tmp_path: Path,
) -> None:
    dataset = load_rabies_method_sensitivity_panel_dataset()
    altered_metadata = tmp_path / "metadata.csv"
    altered_metadata.write_text(
        dataset.metadata_path.read_text(encoding="utf-8")
        + "bad_taxon,bad_host,bad_region\n",
        encoding="utf-8",
    )

    report = audit_rabies_method_sensitivity_workflow_bundle(
        dataset.reference_output_root,
        sequences_path=dataset.sequences_path,
        metadata_path=altered_metadata,
    )

    assert report.all_passed is False
    metadata_check = next(
        row for row in report.checks if row.check_id == "input-checksum:metadata.csv"
    )
    assert metadata_check.status == "failed"


def test_write_rabies_method_sensitivity_reproducibility_artifacts(
    tmp_path: Path,
) -> None:
    dataset = load_rabies_method_sensitivity_panel_dataset()
    report = audit_rabies_method_sensitivity_workflow_bundle(
        dataset.reference_output_root,
        sequences_path=dataset.sequences_path,
        metadata_path=dataset.metadata_path,
    )

    checks_path = write_rabies_method_sensitivity_reproducibility_checks_table(
        tmp_path / "reproducibility-checks.tsv",
        report,
    )
    variants_path = write_rabies_method_sensitivity_variant_audit_table(
        tmp_path / "reproducibility-variants.tsv",
        report,
    )
    audit_path = write_rabies_method_sensitivity_reproducibility_audit_json(
        tmp_path / "reproducibility-audit.json",
        report,
    )

    assert checks_path.read_text(encoding="utf-8").startswith(
        "check_id\tsurface\tstatus\texpected\tobserved\tdetail\n"
    )
    assert variants_path.read_text(encoding="utf-8").startswith(
        "variant_id\tstatus\toutput_file_count\toutput_byte_count\toutput_digest\tmissing_required_files\tunexpected_files\tissues\n"
    )
    payload = json.loads(audit_path.read_text(encoding="utf-8"))
    assert payload["all_passed"] is True
    assert payload["variant_count"] == 4
