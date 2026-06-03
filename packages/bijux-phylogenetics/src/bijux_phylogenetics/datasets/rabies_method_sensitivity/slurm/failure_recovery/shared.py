from __future__ import annotations

import csv
import json
from pathlib import Path

_CONFIG_FILENAME = "workflow-config.resolved.json"
_SLURM_JOB_STATUS_FILENAME = "slurm-job-status.tsv"
_SLURM_PARTITION_STATUS_FILENAME = "slurm-partition-status.tsv"
_SLURM_WORKFLOW_STATUS_FILENAME = "slurm-workflow-status.json"
_TERMINAL_FAILURE_CODES = {
    "parallel_variant_failed": "task_failure",
    "engine_process_timeout": "task_timeout",
    "engine_output_missing": "incomplete_outputs",
    "engine_workflow_running_marker_invalid": "stale_running_marker",
}


def _parse_task_log(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if ": " not in raw_line:
            continue
        key, value = raw_line.split(": ", 1)
        values[key.strip()] = value.strip()
    return values


def _normalize_optional(value: object | None) -> str | None:
    if value is None:
        return None
    text = str(value)
    return None if text == "" else text


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
