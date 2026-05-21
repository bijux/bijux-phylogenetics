from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmTreeRetentionCheckRow:
    """One machine-readable consistency check behind the tree-retention policy."""

    check_id: str
    surface: str
    status: str
    expected: str
    observed: str
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmTreeRetentionFileRow:
    """One tree-bearing file and the safe retention policy derived for it."""

    variant_id: str
    relative_path: str
    artifact_scope: str
    tree_count: int
    byte_count: int
    thinning_policy: str
    thinning_interval: int
    retained_tree_count: int
    compression_policy: str
    recommended_suffix: str
    issue_count: int
    issues: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmTreeRetentionReport:
    """One reviewer-facing tree-retention policy for the governed rabies batch workflow."""

    dataset_id: str
    workflow_prefix: str
    bundle_root: Path
    overall_policy_status: str
    variant_count: int
    file_count: int
    check_count: int
    failed_check_count: int
    tree_artifact_file_count: int
    tree_set_file_count: int
    posterior_sample_file_count: int
    thinning_recommended_file_count: int
    thinning_required_file_count: int
    compression_recommended_file_count: int
    compression_required_file_count: int
    total_tree_count: int
    total_tree_byte_count: int
    largest_tree_set_path: str
    largest_tree_set_tree_count: int
    global_issue_count: int
    global_issues: tuple[str, ...]
    checks: tuple[RabiesMethodSensitivitySlurmTreeRetentionCheckRow, ...]
    files: tuple[RabiesMethodSensitivitySlurmTreeRetentionFileRow, ...]
