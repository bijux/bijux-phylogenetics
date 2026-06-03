from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmStorageAssumptionRow:
    """One explicit assumption behind the preflight storage estimate."""

    assumption_id: str
    parameter: str
    value: str
    rationale: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmStorageCategoryRow:
    """One storage category rolled up across variant-scoped and shared files."""

    category_id: str
    category_label: str
    variant_file_count: int
    workflow_file_count: int
    total_file_count: int
    variant_byte_count: int
    workflow_byte_count: int
    total_byte_count: int
    estimated_storage_mib: int
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmStorageVariantRow:
    """One variant-scoped retained storage estimate broken down by category."""

    variant_id: str
    output_file_count: int
    log_file_count: int
    tree_file_count: int
    posterior_sample_file_count: int
    report_file_count: int
    total_file_count: int
    output_byte_count: int
    log_byte_count: int
    tree_byte_count: int
    posterior_sample_byte_count: int
    report_byte_count: int
    total_byte_count: int
    estimated_storage_mib: int


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmStorageReport:
    """One preflight retained-storage estimate for the governed rabies batch workflow."""

    dataset_id: str
    workflow_prefix: str
    bundle_root: Path
    variant_count: int
    total_file_count: int
    variant_scoped_file_count: int
    workflow_shared_file_count: int
    total_byte_count: int
    total_estimated_storage_mib: int
    variant_scoped_byte_count: int
    workflow_shared_byte_count: int
    output_byte_count: int
    log_byte_count: int
    tree_byte_count: int
    posterior_sample_byte_count: int
    report_byte_count: int
    largest_variant_id: str
    largest_variant_total_byte_count: int
    assumptions: tuple[RabiesMethodSensitivitySlurmStorageAssumptionRow, ...]
    categories: tuple[RabiesMethodSensitivitySlurmStorageCategoryRow, ...]
    variants: tuple[RabiesMethodSensitivitySlurmStorageVariantRow, ...]
