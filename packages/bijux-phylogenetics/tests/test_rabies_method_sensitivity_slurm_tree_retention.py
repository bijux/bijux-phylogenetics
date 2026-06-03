from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm.tree_retention import (
    build_rabies_method_sensitivity_slurm_tree_retention_report,
    write_rabies_method_sensitivity_slurm_tree_retention_checks_table,
    write_rabies_method_sensitivity_slurm_tree_retention_files_table,
    write_rabies_method_sensitivity_slurm_tree_retention_html_report,
    write_rabies_method_sensitivity_slurm_tree_retention_summary_json,
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
            "variants": [
                {"variant_id": "auto-gap-threshold"},
                {"variant_id": "ginsi-gappyout"},
            ],
        },
    )
    _write_tsv(
        bundle_root / "slurm-storage-categories.tsv",
        fieldnames=("category_id", "total_file_count", "total_byte_count"),
        rows=[
            {"category_id": "logs", "total_file_count": 2, "total_byte_count": 4096},
            {"category_id": "outputs", "total_file_count": 4, "total_byte_count": 8192},
            {
                "category_id": "posterior_samples",
                "total_file_count": 1,
                "total_byte_count": 151000,
            },
            {
                "category_id": "reports",
                "total_file_count": 8,
                "total_byte_count": 16384,
            },
            {"category_id": "trees", "total_file_count": 2, "total_byte_count": 152800},
        ],
    )
    _write_json(
        bundle_root / "slurm-output-explosion-report.json",
        {
            "total_tree_byte_count": 152800,
            "total_posterior_sample_byte_count": 151000,
        },
    )
    small_tree = bundle_root / "variants" / "auto-gap-threshold" / "iqtree-support.nwk"
    small_tree.parent.mkdir(parents=True, exist_ok=True)
    small_tree.write_text("(a:0.1,b:0.2);\n", encoding="utf-8")
    big_tree_set = (
        bundle_root / "variants" / "ginsi-gappyout" / "posterior-samples.trees"
    )
    big_tree_set.parent.mkdir(parents=True, exist_ok=True)
    trees = ["(a:0.1,b:0.2);\n" for _ in range(12000)]
    big_tree_set.write_text("".join(trees), encoding="utf-8")


def test_build_rabies_method_sensitivity_slurm_tree_retention_report_flags_large_tree_sets(
    tmp_path: Path,
) -> None:
    _write_synthetic_bundle(tmp_path)

    report = build_rabies_method_sensitivity_slurm_tree_retention_report(tmp_path)

    assert report.dataset_id == "rabies_method_sensitivity_panel"
    assert report.overall_policy_status == "required"
    assert report.file_count == 2
    assert report.tree_set_file_count == 1
    assert report.posterior_sample_file_count == 1
    assert report.thinning_required_file_count == 1
    assert report.compression_required_file_count == 1
    large_row = next(
        row
        for row in report.files
        if row.relative_path.endswith("posterior-samples.trees")
    )
    assert large_row.artifact_scope == "posterior_sample"
    assert large_row.tree_count == 12000
    assert large_row.thinning_policy == "thin_required"
    assert large_row.compression_policy == "compress_required"


def test_build_rabies_method_sensitivity_slurm_tree_retention_report_is_no_action_for_packaged_outputs() -> (
    None
):
    from bijux_phylogenetics.datasets import (
        load_rabies_method_sensitivity_panel_dataset,
    )

    dataset = load_rabies_method_sensitivity_panel_dataset()
    report = build_rabies_method_sensitivity_slurm_tree_retention_report(
        dataset.reference_output_root
    )

    assert report.variant_count == 4
    assert report.overall_policy_status == "no_action"
    assert report.tree_set_file_count == 0
    assert report.thinning_required_file_count == 0
    assert report.compression_required_file_count == 0


def test_write_rabies_method_sensitivity_slurm_tree_retention_artifacts(
    tmp_path: Path,
) -> None:
    _write_synthetic_bundle(tmp_path)
    report = build_rabies_method_sensitivity_slurm_tree_retention_report(tmp_path)

    checks_path = write_rabies_method_sensitivity_slurm_tree_retention_checks_table(
        tmp_path / "slurm-tree-retention-checks.tsv",
        report,
    )
    files_path = write_rabies_method_sensitivity_slurm_tree_retention_files_table(
        tmp_path / "slurm-tree-retention-files.tsv",
        report,
    )
    summary_path = write_rabies_method_sensitivity_slurm_tree_retention_summary_json(
        tmp_path / "slurm-tree-retention-policy.json",
        report,
    )
    html_path = write_rabies_method_sensitivity_slurm_tree_retention_html_report(
        tmp_path / "slurm-tree-retention-policy.html",
        report,
    )

    assert "check_id" in checks_path.read_text(encoding="utf-8")
    assert "thinning_policy" in files_path.read_text(encoding="utf-8")
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["overall_policy_status"] == "required"
    html = html_path.read_text(encoding="utf-8")
    assert "Rabies Slurm Tree Retention Report" in html
    assert "slurm-tree-retention-policy.json" in html
