from __future__ import annotations

import csv
import json
from pathlib import Path

_CONFIG_FILENAME = "workflow-config.resolved.json"
_WORKFLOW_SUMMARY_FILENAME = "workflow-summary.tsv"
_VARIANT_SUMMARY_FILENAME = "variant-summary.tsv"
_PREPROCESSING_COMPARISONS_FILENAME = "preprocessing-rooted-comparisons.tsv"
_STABLE_CLADES_FILENAME = "stable-clades.tsv"
_CHANGED_CLADES_FILENAME = "changed-clades.tsv"
_CONCLUSION_SUMMARY_FILENAME = "method-conclusion-summary.tsv"
_SLURM_JOB_STATUS_FILENAME = "slurm-job-status.tsv"
_SLURM_OUTPUT_FRESHNESS_FILENAME = "slurm-output-freshness.tsv"
_SLURM_JOB_EVIDENCE_FILENAME = "slurm-job-evidence.tsv"


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
