from __future__ import annotations

import csv
import json
from pathlib import Path

_CONFIG_FILENAME = "workflow-config.resolved.json"
_SLURM_JOB_PLAN_FILENAME = "slurm-job-plan.tsv"
_SLURM_STORAGE_CATEGORIES_FILENAME = "slurm-storage-categories.tsv"
_SLURM_STORAGE_VARIANTS_FILENAME = "slurm-storage-variants.tsv"
_SLURM_STORAGE_SUMMARY_FILENAME = "slurm-storage-report.json"
_CATEGORY_IDS = ["logs", "outputs", "posterior_samples", "reports", "trees"]
_MEBIBYTE = 1024 * 1024
_OUTPUT_WARNING_MIB = 128
_OUTPUT_HIGH_MIB = 512
_STORAGE_WARNING_MIB = 256
_STORAGE_HIGH_MIB = 1024
_TREE_WARNING_BYTES = 64 * _MEBIBYTE
_TREE_HIGH_BYTES = 256 * _MEBIBYTE
_TREE_WARNING_FILES = 64
_TREE_HIGH_FILES = 512
_POSTERIOR_WARNING_BYTES = 64 * _MEBIBYTE
_POSTERIOR_HIGH_BYTES = 512 * _MEBIBYTE
_POSTERIOR_WARNING_FILES = 64
_POSTERIOR_HIGH_FILES = 512
_REPORT_WARNING_BYTES = 128 * _MEBIBYTE
_REPORT_HIGH_BYTES = 512 * _MEBIBYTE
_TOTAL_OUTPUT_WARNING_MIB = 512
_TOTAL_OUTPUT_HIGH_MIB = 2048
_TOTAL_STORAGE_WARNING_MIB = 1024
_TOTAL_STORAGE_HIGH_MIB = 4096
_DOMINANT_VARIANT_OUTPUT_SHARE = 0.65
_DOMINANT_VARIANT_OUTPUT_MIB = 128
_BOOTSTRAP_WARNING_REPLICATES = 1000
_BOOTSTRAP_HIGH_REPLICATES = 5000


def _format_float(value: float) -> str:
    return f"{value:.4f}".rstrip("0").rstrip(".")


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


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
            writer.writerow(row)
    return path
