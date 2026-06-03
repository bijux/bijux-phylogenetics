from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .shared import (
    _CONFIG_FILENAME,
    _SLURM_JOB_STATUS_FILENAME,
    _SLURM_PARTITION_STATUS_FILENAME,
    _SLURM_WORKFLOW_STATUS_FILENAME,
    _load_json,
    _read_tsv_rows,
)


@dataclass(frozen=True, slots=True)
class FailureRecoveryInputs:
    bundle_root: Path
    config: dict[str, object]
    configured_variant_ids: list[str]
    job_status_rows: list[dict[str, str]]
    partition_status_rows: list[dict[str, str]]
    workflow_status: dict[str, object]
    checks: tuple[tuple[str, str, bool, object, object, str], ...]


def load_failure_recovery_inputs(bundle_root: Path) -> FailureRecoveryInputs:
    bundle_root = bundle_root.resolve()
    config = _load_json(bundle_root / _CONFIG_FILENAME)
    job_status_rows = _read_tsv_rows(bundle_root / _SLURM_JOB_STATUS_FILENAME)
    partition_status_rows = _read_tsv_rows(
        bundle_root / _SLURM_PARTITION_STATUS_FILENAME
    )
    workflow_status = _load_json(bundle_root / _SLURM_WORKFLOW_STATUS_FILENAME)

    configured_variant_ids = sorted(
        str(row["variant_id"]) for row in list(config.get("variants", []))
    )
    observed_variant_ids = sorted(str(row["variant_id"]) for row in job_status_rows)
    checks = (
        (
            "job-status:variant-coverage",
            "job-status",
            observed_variant_ids == configured_variant_ids,
            configured_variant_ids,
            observed_variant_ids,
            "job-status rows cover the configured variant ids",
        ),
    )
    return FailureRecoveryInputs(
        bundle_root=bundle_root,
        config=config,
        configured_variant_ids=configured_variant_ids,
        job_status_rows=job_status_rows,
        partition_status_rows=partition_status_rows,
        workflow_status=workflow_status,
        checks=checks,
    )
