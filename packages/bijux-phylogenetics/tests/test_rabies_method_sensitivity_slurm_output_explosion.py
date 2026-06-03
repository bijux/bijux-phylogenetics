from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm.output_explosion import (
    build_rabies_method_sensitivity_slurm_output_explosion_report,
    write_rabies_method_sensitivity_slurm_output_explosion_checks_table,
    write_rabies_method_sensitivity_slurm_output_explosion_html_report,
    write_rabies_method_sensitivity_slurm_output_explosion_summary_json,
    write_rabies_method_sensitivity_slurm_output_explosion_variants_table,
)


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _write_tsv(
    path: Path,
    *,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_synthetic_bundle(bundle_root: Path) -> None:
    _write_json(
        bundle_root / "workflow-config.resolved.json",
        {
            "dataset_id": "rabies_method_sensitivity_panel",
            "workflow_prefix": "rabies-method-sensitivity-panel",
            "bootstrap_replicates": 1000,
            "variants": [
                {"variant_id": "auto-gap-threshold"},
                {"variant_id": "ginsi-gappyout"},
            ],
        },
    )
    _write_tsv(
        bundle_root / "slurm-job-plan.tsv",
        fieldnames=(
            "variant_id",
            "estimated_output_mib",
        ),
        rows=[
            {"variant_id": "auto-gap-threshold", "estimated_output_mib": 16},
            {"variant_id": "ginsi-gappyout", "estimated_output_mib": 640},
        ],
    )
    _write_tsv(
        bundle_root / "slurm-storage-categories.tsv",
        fieldnames=(
            "category_id",
            "total_file_count",
            "total_byte_count",
        ),
        rows=[
            {"category_id": "logs", "total_file_count": 2, "total_byte_count": 4096},
            {
                "category_id": "outputs",
                "total_file_count": 4,
                "total_byte_count": 2_097_152,
            },
            {
                "category_id": "posterior_samples",
                "total_file_count": 128,
                "total_byte_count": 134_217_728,
            },
            {
                "category_id": "reports",
                "total_file_count": 12,
                "total_byte_count": 16_777_216,
            },
            {
                "category_id": "trees",
                "total_file_count": 96,
                "total_byte_count": 83_886_080,
            },
        ],
    )
    _write_tsv(
        bundle_root / "slurm-storage-variants.tsv",
        fieldnames=(
            "variant_id",
            "tree_file_count",
            "tree_byte_count",
            "posterior_sample_file_count",
            "posterior_sample_byte_count",
            "report_byte_count",
            "estimated_storage_mib",
            "total_byte_count",
        ),
        rows=[
            {
                "variant_id": "auto-gap-threshold",
                "tree_file_count": 4,
                "tree_byte_count": 2048,
                "posterior_sample_file_count": 0,
                "posterior_sample_byte_count": 0,
                "report_byte_count": 4096,
                "estimated_storage_mib": 1,
                "total_byte_count": 8192,
            },
            {
                "variant_id": "ginsi-gappyout",
                "tree_file_count": 92,
                "tree_byte_count": 83_884_032,
                "posterior_sample_file_count": 128,
                "posterior_sample_byte_count": 134_217_728,
                "report_byte_count": 16_773_120,
                "estimated_storage_mib": 300,
                "total_byte_count": 234_881_024,
            },
        ],
    )
    _write_json(
        bundle_root / "slurm-storage-report.json",
        {
            "variant_count": 2,
            "total_byte_count": 236_978_176,
            "total_estimated_storage_mib": 227,
            "largest_variant_id": "ginsi-gappyout",
            "largest_variant_total_byte_count": 234_881_024,
        },
    )


def test_build_rabies_method_sensitivity_slurm_output_explosion_report_detects_high_risk_variant(
    tmp_path: Path,
) -> None:
    _write_synthetic_bundle(tmp_path)

    report = build_rabies_method_sensitivity_slurm_output_explosion_report(tmp_path)

    assert report.dataset_id == "rabies_method_sensitivity_panel"
    assert report.overall_risk_status == "high"
    assert report.warning_variant_count == 0
    assert report.high_risk_variant_count == 1
    assert report.total_estimated_output_mib == 656
    assert report.total_posterior_sample_byte_count == 134_217_728
    risky_variant = next(
        row for row in report.variants if row.variant_id == "ginsi-gappyout"
    )
    assert risky_variant.risk_status == "high"
    assert risky_variant.issue_count > 0
    assert any("posterior sample" in issue for issue in risky_variant.issues)


def test_build_rabies_method_sensitivity_slurm_output_explosion_report_is_low_risk_for_packaged_outputs() -> (
    None
):
    from bijux_phylogenetics.datasets import (
        load_rabies_method_sensitivity_panel_dataset,
    )

    dataset = load_rabies_method_sensitivity_panel_dataset()
    report = build_rabies_method_sensitivity_slurm_output_explosion_report(
        dataset.reference_output_root
    )

    assert report.variant_count == 4
    assert report.overall_risk_status == "low"
    assert report.high_risk_variant_count == 0
    assert report.total_posterior_sample_byte_count == 0


def test_write_rabies_method_sensitivity_slurm_output_explosion_artifacts(
    tmp_path: Path,
) -> None:
    _write_synthetic_bundle(tmp_path)
    report = build_rabies_method_sensitivity_slurm_output_explosion_report(tmp_path)

    checks_path = write_rabies_method_sensitivity_slurm_output_explosion_checks_table(
        tmp_path / "slurm-output-explosion-checks.tsv",
        report,
    )
    variants_path = (
        write_rabies_method_sensitivity_slurm_output_explosion_variants_table(
            tmp_path / "slurm-output-explosion-variants.tsv",
            report,
        )
    )
    summary_path = write_rabies_method_sensitivity_slurm_output_explosion_summary_json(
        tmp_path / "slurm-output-explosion-report.json",
        report,
    )
    html_path = write_rabies_method_sensitivity_slurm_output_explosion_html_report(
        tmp_path / "slurm-output-explosion-report.html",
        report,
    )

    assert "check_id" in checks_path.read_text(encoding="utf-8")
    assert "risk_status" in variants_path.read_text(encoding="utf-8")
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["overall_risk_status"] == "high"
    html = html_path.read_text(encoding="utf-8")
    assert "Rabies Slurm Output Explosion Report" in html
    assert "slurm-output-explosion-report.json" in html
