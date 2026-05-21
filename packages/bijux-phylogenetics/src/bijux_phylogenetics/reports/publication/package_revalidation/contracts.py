from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class PublicationPackageRevalidationArtifactRow:
    """One reviewer-facing revalidation row for a declared package artifact."""

    artifact_scope: str
    section: str
    kind: str
    relative_path: str
    status: str
    expected_sha256: str | None
    observed_sha256: str | None
    expected_size_bytes: int | None
    observed_size_bytes: int | None
    detail: str


@dataclass(frozen=True, slots=True)
class PublicationPackageRevalidationCheckRow:
    """One reviewer-facing package revalidation decision row."""

    section: str
    check_id: str
    status: str
    summary: str
    evidence: str
    artifact_path: str


@dataclass(slots=True)
class PublicationPackageRevalidationResult:
    """Written revalidation artifacts for one stored publication package."""

    output_root: Path
    manifest_path: Path
    package_root: Path
    report_kind: str
    artifact_table_path: Path
    check_table_path: Path
    summary_path: Path
    report_path: Path
    artifact_rows: list[PublicationPackageRevalidationArtifactRow]
    check_rows: list[PublicationPackageRevalidationCheckRow]
    matched_artifact_count: int
    missing_artifact_count: int
    checksum_mismatch_count: int
    size_mismatch_count: int
    unexpected_file_count: int
    blocked_check_count: int
    risk_check_count: int
    all_original_artifacts_match: bool
    overall_revalidation_status: str
