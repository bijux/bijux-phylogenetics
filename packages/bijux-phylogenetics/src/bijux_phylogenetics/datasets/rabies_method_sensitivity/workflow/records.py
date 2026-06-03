from __future__ import annotations

import json
from pathlib import Path

from ..models import (
    RabiesMethodSensitivityPanelDataset,
    RabiesMethodSensitivityTaskRecord,
    RabiesMethodSensitivityVariant,
)
from ..shared import _WORKFLOW_PREFIX


def _workflow_execution_record_path(output_root: Path) -> Path:
    """Return the workflow-root execution record for one packaged sensitivity run."""
    return output_root / f"{_WORKFLOW_PREFIX}.run.json"


def _write_workflow_execution_record(
    path: Path,
    *,
    dataset: RabiesMethodSensitivityPanelDataset,
    execution_mode: str,
    parallel_workers: int,
    task_records: tuple[RabiesMethodSensitivityTaskRecord, ...],
    started_at_utc: str,
    ended_at_utc: str,
    status: str,
) -> Path:
    """Persist one execution record for one packaged sensitivity workflow run."""
    payload = {
        "dataset_id": dataset.dataset_id,
        "workflow": "rabies_method_sensitivity_panel",
        "workflow_prefix": dataset.workflow_prefix,
        "status": status,
        "started_at_utc": started_at_utc,
        "ended_at_utc": ended_at_utc,
        "parallel_workers": parallel_workers,
        "execution_mode": execution_mode,
        "variant_count": len(dataset.variants),
        "selected_variant_ids": [variant.variant_id for variant in dataset.variants],
        "successful_variants": [
            record.variant_id for record in task_records if record.status == "succeeded"
        ],
        "failed_variants": [
            record.variant_id for record in task_records if record.status != "succeeded"
        ],
        "task_records": [
            {
                "variant_id": record.variant_id,
                "label": record.label,
                "status": record.status,
                "execution_mode": record.execution_mode,
                "log_path": Path("parallel-logs", record.log_path.name).as_posix(),
                "output_root": Path("variants", record.variant_id).as_posix(),
                "error_code": record.error_code,
                "error_message": record.error_message,
            }
            for record in task_records
        ],
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _write_task_log(
    path: Path,
    *,
    variant: RabiesMethodSensitivityVariant,
    execution_mode: str,
    status: str,
    output_root: Path,
    error_code: str | None,
    error_message: str | None,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"variant_id: {variant.variant_id}",
        f"label: {variant.label}",
        f"execution_mode: {execution_mode}",
        f"alignment_mode: {variant.alignment_mode}",
        f"trimming_mode: {variant.trimming_mode}",
        f"trim_gap_threshold: {_format_float(variant.trim_gap_threshold)}",
        f"status: {status}",
        f"output_root: {output_root.as_posix()}",
    ]
    if error_code is not None:
        lines.append(f"error_code: {error_code}")
    if error_message is not None:
        lines.append(f"error_message: {error_message}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _format_float(value: float) -> str:
    return format(value, ".12g")
