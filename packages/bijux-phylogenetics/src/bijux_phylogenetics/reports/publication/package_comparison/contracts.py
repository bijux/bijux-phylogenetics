from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class PublicationPackageComparisonArtifactRow:
    """One artifact-level difference row between two stored study packages."""

    section: str
    kind: str
    relative_path: str
    status: str
    left_sha256: str | None
    right_sha256: str | None
    left_size_bytes: int | None
    right_size_bytes: int | None
    detail: str


@dataclass(frozen=True, slots=True)
class PublicationPackageComparisonCheckRow:
    """One reviewer-facing comparison decision over two package versions."""

    section: str
    check_id: str
    status: str
    summary: str
    evidence: str
    left_artifact_path: str
    right_artifact_path: str


@dataclass(slots=True)
class PublicationPackageComparisonResult:
    """Written comparison artifacts for two stored publication packages."""

    output_root: Path
    left_manifest_path: Path
    right_manifest_path: Path
    left_package_root: Path
    right_package_root: Path
    report_kind: str
    dataset_id: str
    artifact_table_path: Path
    check_table_path: Path
    summary_path: Path
    report_path: Path
    artifact_rows: list[PublicationPackageComparisonArtifactRow]
    check_rows: list[PublicationPackageComparisonCheckRow]
    same_artifact_count: int
    changed_artifact_count: int
    left_only_artifact_count: int
    right_only_artifact_count: int
    config_difference_count: int
    sequence_left_only_count: int
    sequence_right_only_count: int
    accession_left_only_count: int
    accession_right_only_count: int
    alignment_difference_count: int
    figure_or_report_difference_count: int
    scientific_finding_difference_count: int
    overall_comparison_status: str
