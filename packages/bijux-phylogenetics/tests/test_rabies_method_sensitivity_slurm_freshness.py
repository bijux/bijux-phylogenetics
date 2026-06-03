from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path

from bijux_phylogenetics.datasets.rabies_method_sensitivity import (
    load_rabies_method_sensitivity_panel_dataset,
)
from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm.freshness import (
    build_rabies_method_sensitivity_slurm_output_freshness_report,
    write_rabies_method_sensitivity_slurm_output_freshness_checks_table,
    write_rabies_method_sensitivity_slurm_output_freshness_json,
    write_rabies_method_sensitivity_slurm_output_freshness_table,
)


def test_build_rabies_method_sensitivity_slurm_output_freshness_report_passes_on_packaged_outputs() -> (
    None
):
    dataset = load_rabies_method_sensitivity_panel_dataset()

    report = build_rabies_method_sensitivity_slurm_output_freshness_report(
        dataset.reference_output_root,
        dataset=dataset,
    )

    assert report.all_outputs_fresh is True
    assert report.failed_check_count == 0
    assert report.fresh_job_count == 4
    assert report.stale_job_count == 0
    assert all(row.freshness_status == "fresh" for row in report.jobs)


def test_build_rabies_method_sensitivity_slurm_output_freshness_report_detects_input_drift(
    tmp_path: Path,
) -> None:
    dataset = load_rabies_method_sensitivity_panel_dataset()
    drifted_metadata = tmp_path / "metadata.csv"
    drifted_metadata.write_text(
        dataset.metadata_path.read_text(encoding="utf-8")
        + "drift_taxon,drift_host,drift_region\n",
        encoding="utf-8",
    )
    drifted_dataset = replace(dataset, metadata_path=drifted_metadata)

    report = build_rabies_method_sensitivity_slurm_output_freshness_report(
        dataset.reference_output_root,
        dataset=drifted_dataset,
    )

    assert report.all_outputs_fresh is False
    assert report.stale_job_count == 4
    metadata_check = next(
        row for row in report.checks if row.check_id == "input-checksum:metadata.csv"
    )
    assert metadata_check.status == "failed"
    assert all(
        "input-checksum:metadata.csv" in row.stale_reason_codes for row in report.jobs
    )


def test_build_rabies_method_sensitivity_slurm_output_freshness_report_detects_variant_setting_drift(
    tmp_path: Path,
) -> None:
    dataset = load_rabies_method_sensitivity_panel_dataset()
    changed_variant = replace(
        dataset.variants[0],
        trim_gap_threshold=dataset.variants[0].trim_gap_threshold + 0.05,
    )
    drifted_dataset = replace(
        dataset,
        variants=(changed_variant, *dataset.variants[1:]),
    )

    report = build_rabies_method_sensitivity_slurm_output_freshness_report(
        dataset.reference_output_root,
        dataset=drifted_dataset,
    )

    assert report.all_outputs_fresh is False
    stale_rows = [row for row in report.jobs if row.freshness_status == "stale"]
    assert [row.variant_id for row in stale_rows] == [dataset.variants[0].variant_id]
    assert stale_rows[0].variant_settings_match is False
    assert (
        f"variant-setting:{dataset.variants[0].variant_id}"
        in stale_rows[0].stale_reason_codes
    )


def test_write_rabies_method_sensitivity_slurm_output_freshness_artifacts(
    tmp_path: Path,
) -> None:
    dataset = load_rabies_method_sensitivity_panel_dataset()
    report = build_rabies_method_sensitivity_slurm_output_freshness_report(
        dataset.reference_output_root,
        dataset=dataset,
    )

    table_path = write_rabies_method_sensitivity_slurm_output_freshness_table(
        tmp_path / "slurm-output-freshness.tsv",
        report,
    )
    checks_path = write_rabies_method_sensitivity_slurm_output_freshness_checks_table(
        tmp_path / "slurm-output-freshness-checks.tsv",
        report,
    )
    summary_path = write_rabies_method_sensitivity_slurm_output_freshness_json(
        tmp_path / "slurm-output-freshness.json",
        report,
    )

    assert "freshness_status" in table_path.read_text(encoding="utf-8")
    assert "check_id" in checks_path.read_text(encoding="utf-8")
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["bundle_root"] == "."
    assert payload["fresh_job_count"] == 4
    assert payload["jobs"][0]["freshness_status"] == "fresh"
