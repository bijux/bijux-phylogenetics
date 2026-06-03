from __future__ import annotations

from pathlib import Path

from ..support import artifact_kind, checksum
from .contracts import PublicationPackageRevalidationArtifactRow


def artifact_row(
    *,
    artifact_scope: str,
    section: str,
    relative_path: str,
    expected_sha256: str | None,
    expected_size_bytes: int | None,
    observed_path: Path,
) -> PublicationPackageRevalidationArtifactRow:
    if not observed_path.exists():
        return PublicationPackageRevalidationArtifactRow(
            artifact_scope=artifact_scope,
            section=section,
            kind=artifact_kind(relative_path),
            relative_path=relative_path,
            status="blocked",
            expected_sha256=expected_sha256,
            observed_sha256=None,
            expected_size_bytes=expected_size_bytes,
            observed_size_bytes=None,
            detail="declared package artifact is missing",
        )
    observed_size = observed_path.stat().st_size
    observed_sha256 = checksum(observed_path)
    mismatches: list[str] = []
    if expected_size_bytes is not None and observed_size != expected_size_bytes:
        mismatches.append(
            f"stored size {expected_size_bytes} bytes != observed size {observed_size} bytes"
        )
    if expected_sha256 is not None and observed_sha256 != expected_sha256:
        mismatches.append("stored checksum does not match observed checksum")
    return PublicationPackageRevalidationArtifactRow(
        artifact_scope=artifact_scope,
        section=section,
        kind=artifact_kind(relative_path),
        relative_path=relative_path,
        status="blocked" if mismatches else "pass",
        expected_sha256=expected_sha256,
        observed_sha256=observed_sha256,
        expected_size_bytes=expected_size_bytes,
        observed_size_bytes=observed_size,
        detail="; ".join(mismatches)
        if mismatches
        else "artifact matches stored package record",
    )


def unexpected_files(
    *,
    package_root: Path,
    expected_relative_paths: set[str],
    output_root: Path,
    ignored_prefixes: tuple[str, ...],
) -> list[str]:
    unexpected: list[str] = []
    output_root_resolved = output_root.resolve()
    for path in sorted(package_root.rglob("*")):
        if not path.is_file():
            continue
        if output_root_resolved in path.resolve().parents:
            continue
        relative_path = path.relative_to(package_root).as_posix()
        if any(relative_path.startswith(prefix) for prefix in ignored_prefixes):
            continue
        if relative_path not in expected_relative_paths:
            unexpected.append(relative_path)
    return unexpected
