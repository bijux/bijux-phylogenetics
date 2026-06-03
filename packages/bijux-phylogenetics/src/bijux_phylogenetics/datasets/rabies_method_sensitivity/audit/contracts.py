from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityReproducibilityCheckRow:
    """One machine-readable pass/fail check within the bundle audit."""

    check_id: str
    surface: str
    status: str
    expected: str
    observed: str
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityVariantAuditRow:
    """One per-variant provenance and file-inventory summary."""

    variant_id: str
    status: str
    output_file_count: int
    output_byte_count: int
    output_digest: str
    missing_required_files: tuple[str, ...]
    unexpected_files: tuple[str, ...]
    issues: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityReproducibilityAuditReport:
    """One reviewer-facing reproducibility audit for the workflow bundle."""

    dataset_id: str
    bundle_root: Path
    workflow_manifest_path: Path
    report_manifest_path: Path
    config_path: Path
    sequences_path: Path
    metadata_path: Path
    all_passed: bool
    check_count: int
    failed_check_count: int
    variant_count: int
    failed_variant_count: int
    checks: tuple[RabiesMethodSensitivityReproducibilityCheckRow, ...]
    variants: tuple[RabiesMethodSensitivityVariantAuditRow, ...]
