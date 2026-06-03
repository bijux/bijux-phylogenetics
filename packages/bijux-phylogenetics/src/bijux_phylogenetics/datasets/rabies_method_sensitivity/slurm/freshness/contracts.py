from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityOutputFreshnessCheckRow:
    """One freshness check over bundle inputs or output-affecting settings."""

    check_id: str
    surface: str
    scope: str
    status: str
    expected: str
    observed: str
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmOutputFreshnessRow:
    """One per-job output-freshness classification against current inputs/settings."""

    partition_id: str
    array_index: int
    variant_id: str
    freshness_status: str
    inputs_match: bool
    workflow_settings_match: bool
    variant_settings_match: bool
    stale_reason_count: int
    stale_reason_codes: tuple[str, ...]
    stale_reason_detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmOutputFreshnessReport:
    """One reviewer-facing freshness summary for rabies batch outputs."""

    dataset_id: str
    workflow_prefix: str
    bundle_root: Path
    all_outputs_fresh: bool
    selected_variant_ids: tuple[str, ...]
    check_count: int
    failed_check_count: int
    fresh_job_count: int
    stale_job_count: int
    checks: tuple[RabiesMethodSensitivityOutputFreshnessCheckRow, ...]
    jobs: tuple[RabiesMethodSensitivitySlurmOutputFreshnessRow, ...]
