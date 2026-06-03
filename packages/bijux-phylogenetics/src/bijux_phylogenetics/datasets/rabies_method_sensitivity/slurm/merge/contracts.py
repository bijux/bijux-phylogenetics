from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmMergeCheckRow:
    """One machine-readable merge-consistency check over batch outputs."""

    check_id: str
    surface: str
    status: str
    expected: str
    observed: str
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmMergeVariantRow:
    """One per-variant decision about whether the job can join the merged result."""

    variant_id: str
    merge_status: str
    job_status: str
    output_freshness_status: str
    evidence_status: str
    included_in_merge: bool
    selected_model: str
    serious_conflict_count: int
    rooted_engine_rf_distance: int
    rooted_engine_same_taxa_different_rooting: bool
    issue_count: int
    issues: tuple[str, ...]
    evidence_json_path: str
    evidence_html_path: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmMergeReport:
    """One explicit global merge contract over distributed rabies HPC outputs."""

    dataset_id: str
    workflow_prefix: str
    bundle_root: Path
    merge_status: str
    merge_ready: bool
    expected_variant_count: int
    merged_variant_count: int
    mergeable_variant_count: int
    failed_variant_count: int
    failed_check_count: int
    check_count: int
    stable_clade_count: int
    changed_clade_count: int
    preprocessing_comparison_count: int
    conclusion_count: int
    serious_conflict_variant_count: int
    rooted_engine_change_variant_count: int
    maximum_serious_conflict_count: int
    selected_models: tuple[str, ...]
    checks: tuple[RabiesMethodSensitivitySlurmMergeCheckRow, ...]
    variants: tuple[RabiesMethodSensitivitySlurmMergeVariantRow, ...]
