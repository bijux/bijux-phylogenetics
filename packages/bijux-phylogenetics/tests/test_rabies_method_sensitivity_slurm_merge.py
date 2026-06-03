from __future__ import annotations

import json
from pathlib import Path
import shutil

from bijux_phylogenetics.datasets import (
    build_rabies_method_sensitivity_slurm_merge_report,
    load_rabies_method_sensitivity_panel_dataset,
    write_rabies_method_sensitivity_slurm_merge_checks_table,
    write_rabies_method_sensitivity_slurm_merge_html_report,
    write_rabies_method_sensitivity_slurm_merge_summary_json,
    write_rabies_method_sensitivity_slurm_merge_variants_table,
)


def test_build_rabies_method_sensitivity_slurm_merge_report_passes_on_packaged_outputs() -> (
    None
):
    dataset = load_rabies_method_sensitivity_panel_dataset()

    report = build_rabies_method_sensitivity_slurm_merge_report(
        dataset.reference_output_root
    )

    assert report.merge_ready is True
    assert report.merge_status == "merge-ready"
    assert report.expected_variant_count == 4
    assert report.merged_variant_count == 4
    assert report.failed_variant_count == 0
    assert report.failed_check_count == 0
    assert report.stable_clade_count == 2
    assert report.changed_clade_count == 8
    assert report.maximum_serious_conflict_count == 8
    assert report.selected_models == ("TPM2u+F+G4",)


def test_build_rabies_method_sensitivity_slurm_merge_report_detects_blocked_job(
    tmp_path: Path,
) -> None:
    dataset = load_rabies_method_sensitivity_panel_dataset()
    bundle_root = tmp_path / "workflow"
    shutil.copytree(dataset.reference_output_root, bundle_root)
    job_status_path = bundle_root / "slurm-job-status.tsv"
    text = job_status_path.read_text(encoding="utf-8")
    job_status_path.write_text(
        text.replace("\tcompleted\t", "\tpending\t", 1),
        encoding="utf-8",
    )

    report = build_rabies_method_sensitivity_slurm_merge_report(bundle_root)

    assert report.merge_ready is False
    assert report.merge_status == "merge-blocked"
    assert report.failed_variant_count == 1
    blocked_variant = next(row for row in report.variants if not row.included_in_merge)
    assert blocked_variant.job_status == "pending"
    assert "job-status does not mark the variant as completed" in blocked_variant.issues


def test_write_rabies_method_sensitivity_slurm_merge_artifacts(
    tmp_path: Path,
) -> None:
    dataset = load_rabies_method_sensitivity_panel_dataset()
    report = build_rabies_method_sensitivity_slurm_merge_report(
        dataset.reference_output_root
    )

    checks_path = write_rabies_method_sensitivity_slurm_merge_checks_table(
        tmp_path / "slurm-merge-checks.tsv",
        report,
    )
    variants_path = write_rabies_method_sensitivity_slurm_merge_variants_table(
        tmp_path / "slurm-merge-variants.tsv",
        report,
    )
    summary_path = write_rabies_method_sensitivity_slurm_merge_summary_json(
        tmp_path / "slurm-merge-report.json",
        report,
    )
    html_path = write_rabies_method_sensitivity_slurm_merge_html_report(
        tmp_path / "slurm-merge-report.html",
        report,
    )

    assert checks_path.read_text(encoding="utf-8").startswith(
        "check_id\tsurface\tstatus\texpected\tobserved\tdetail\n"
    )
    assert variants_path.read_text(encoding="utf-8").startswith(
        "variant_id\tmerge_status\tjob_status\toutput_freshness_status\t"
    )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["merge_ready"] is True
    assert summary["merged_variant_count"] == 4
    html = html_path.read_text(encoding="utf-8")
    assert 'href="slurm-merge-checks.tsv"' in html
    assert 'href="slurm-merge-report.json"' in html
