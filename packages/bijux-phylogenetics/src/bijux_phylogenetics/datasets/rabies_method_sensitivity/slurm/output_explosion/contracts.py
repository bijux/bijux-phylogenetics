from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmOutputExplosionCheckRow:
    """One machine-readable consistency check behind the explosion-risk report."""

    check_id: str
    surface: str
    status: str
    expected: str
    observed: str
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmOutputExplosionVariantRow:
    """One per-variant retained-output explosion assessment."""

    variant_id: str
    risk_status: str
    estimated_output_mib: int
    estimated_storage_mib: int
    tree_file_count: int
    tree_byte_count: int
    posterior_sample_file_count: int
    posterior_sample_byte_count: int
    report_byte_count: int
    output_share: float
    issue_count: int
    issues: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmOutputExplosionReport:
    """One retained-output explosion-risk summary over the governed rabies batch workflow."""

    dataset_id: str
    workflow_prefix: str
    bundle_root: Path
    bootstrap_replicates: int
    overall_risk_status: str
    variant_count: int
    check_count: int
    failed_check_count: int
    low_risk_variant_count: int
    warning_variant_count: int
    high_risk_variant_count: int
    global_issue_count: int
    total_estimated_output_mib: int
    total_estimated_storage_mib: int
    total_tree_byte_count: int
    total_tree_file_count: int
    total_posterior_sample_byte_count: int
    total_posterior_sample_file_count: int
    total_report_byte_count: int
    largest_variant_id: str
    largest_variant_output_share: float
    global_issues: tuple[str, ...]
    checks: tuple[RabiesMethodSensitivitySlurmOutputExplosionCheckRow, ...]
    variants: tuple[RabiesMethodSensitivitySlurmOutputExplosionVariantRow, ...]
