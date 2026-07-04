from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm.storage import (
    build_rabies_method_sensitivity_slurm_storage_report,
    write_rabies_method_sensitivity_slurm_storage_categories_table,
    write_rabies_method_sensitivity_slurm_storage_html_report,
    write_rabies_method_sensitivity_slurm_storage_summary_json,
    write_rabies_method_sensitivity_slurm_storage_variants_table,
)


def _write_bytes(path: Path, *, byte_count: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"x" * byte_count)


def _write_config(bundle_root: Path) -> None:
    payload = {
        "dataset_id": "rabies_method_sensitivity_panel",
        "workflow_prefix": "rabies-method-sensitivity-panel",
        "variants": [
            {"variant_id": "auto-gap-threshold"},
            {"variant_id": "ginsi-gappyout"},
        ],
    }
    (bundle_root / "workflow-config.resolved.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )


def test_build_rabies_method_sensitivity_slurm_storage_report_classifies_categories(
    tmp_path: Path,
) -> None:
    _write_config(tmp_path)
    _write_bytes(
        tmp_path / "variants" / "auto-gap-threshold" / "auto-gap-threshold.aln",
        byte_count=120,
    )
    _write_bytes(
        tmp_path / "variants" / "auto-gap-threshold" / "rooted-fasttree.nwk",
        byte_count=40,
    )
    _write_bytes(
        tmp_path / "parallel-logs" / "auto-gap-threshold.log",
        byte_count=55,
    )
    _write_bytes(
        tmp_path / "slurm-job-evidence" / "auto-gap-threshold" / "task.log",
        byte_count=60,
    )
    _write_bytes(
        tmp_path / "variants" / "ginsi-gappyout" / "ginsi-gappyout.aln",
        byte_count=180,
    )
    _write_bytes(
        tmp_path / "variants" / "ginsi-gappyout" / "posterior.trees",
        byte_count=75,
    )
    _write_bytes(
        tmp_path / "rabies-method-sensitivity-report.html",
        byte_count=90,
    )

    report = build_rabies_method_sensitivity_slurm_storage_report(tmp_path)

    assert report.dataset_id == "rabies_method_sensitivity_panel"
    assert report.variant_count == 2
    assert report.total_file_count == 8
    assert report.output_byte_count == 300
    assert report.tree_byte_count == 40
    assert report.log_byte_count == 55
    assert report.posterior_sample_byte_count == 75
    assert report.report_byte_count > 150
    assert report.workflow_shared_byte_count > 90
    categories = {row.category_id: row for row in report.categories}
    assert categories["outputs"].variant_file_count == 2
    assert categories["trees"].variant_byte_count == 40
    assert categories["logs"].variant_file_count == 1
    assert categories["posterior_samples"].total_file_count == 1
    assert categories["reports"].workflow_file_count == 2
    auto_row, ginsi_row = report.variants
    assert report.largest_variant_id == auto_row.variant_id
    assert report.largest_variant_total_byte_count == auto_row.total_byte_count
    assert auto_row.variant_id == "auto-gap-threshold"
    assert auto_row.output_file_count == 1
    assert auto_row.tree_file_count == 1
    assert auto_row.log_file_count == 1
    assert auto_row.report_file_count == 1
    assert ginsi_row.posterior_sample_file_count == 1
    assert ginsi_row.report_file_count == 0


def test_build_rabies_method_sensitivity_slurm_storage_report_keeps_zero_posterior_category_for_packaged_outputs() -> (
    None
):
    from bijux_phylogenetics.datasets import (
        load_rabies_method_sensitivity_panel_dataset,
    )

    dataset = load_rabies_method_sensitivity_panel_dataset()
    report = build_rabies_method_sensitivity_slurm_storage_report(
        dataset.reference_output_root
    )

    assert report.variant_count == 4
    assert report.posterior_sample_byte_count == 0
    assert any(row.category_id == "posterior_samples" for row in report.categories)
    assert report.total_estimated_storage_mib > 0


def test_write_rabies_method_sensitivity_slurm_storage_artifacts(
    tmp_path: Path,
) -> None:
    _write_config(tmp_path)
    _write_bytes(
        tmp_path / "variants" / "auto-gap-threshold" / "auto-gap-threshold.aln",
        byte_count=120,
    )
    _write_bytes(
        tmp_path / "parallel-logs" / "auto-gap-threshold.log",
        byte_count=55,
    )
    _write_bytes(
        tmp_path / "rabies-method-sensitivity-report.html",
        byte_count=90,
    )

    report = build_rabies_method_sensitivity_slurm_storage_report(tmp_path)
    categories_path = write_rabies_method_sensitivity_slurm_storage_categories_table(
        tmp_path / "slurm-storage-categories.tsv",
        report,
    )
    variants_path = write_rabies_method_sensitivity_slurm_storage_variants_table(
        tmp_path / "slurm-storage-variants.tsv",
        report,
    )
    summary_path = write_rabies_method_sensitivity_slurm_storage_summary_json(
        tmp_path / "slurm-storage-report.json",
        report,
    )
    html_path = write_rabies_method_sensitivity_slurm_storage_html_report(
        tmp_path / "slurm-storage-report.html",
        report,
    )

    assert "category_id" in categories_path.read_text(encoding="utf-8")
    assert "variant_id" in variants_path.read_text(encoding="utf-8")
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["dataset_id"] == "rabies_method_sensitivity_panel"
    assert payload["total_estimated_storage_mib"] > 0
    html = html_path.read_text(encoding="utf-8")
    assert "Rabies Slurm Storage Report" in html
    assert "slurm-storage-categories.tsv" in html
