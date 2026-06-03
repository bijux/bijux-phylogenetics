from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.io.fasta import load_fasta_records

from ..support import (
    artifact_kind,
    checksum,
    entry_path,
    mapping,
    read_tsv_rows,
)
from .contracts import PublicationPackageComparisonArtifactRow


def inventory_index(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Index artifact inventory rows by their relative path."""
    return {
        row["relative_path"]: row
        for row in rows
        if row.get("relative_path", "").strip()
    }


def artifact_row(
    relative_path: str,
    left_row: dict[str, str] | None,
    right_row: dict[str, str] | None,
) -> PublicationPackageComparisonArtifactRow:
    """Build one artifact-difference row for one relative path."""
    left_sha = None if left_row is None else left_row.get("sha256", "").strip() or None
    right_sha = (
        None if right_row is None else right_row.get("sha256", "").strip() or None
    )
    left_size = (
        None
        if left_row is None or not left_row.get("size_bytes", "").strip()
        else int(left_row["size_bytes"])
    )
    right_size = (
        None
        if right_row is None or not right_row.get("size_bytes", "").strip()
        else int(right_row["size_bytes"])
    )
    section = (
        left_row.get("section", "").strip()
        if left_row is not None
        else right_row.get("section", "").strip()
        if right_row is not None
        else "artifact"
    )
    kind = (
        left_row.get("kind", "").strip()
        if left_row is not None
        else right_row.get("kind", "").strip()
        if right_row is not None
        else artifact_kind(relative_path)
    )
    if left_row is None:
        return PublicationPackageComparisonArtifactRow(
            section=section,
            kind=kind,
            relative_path=relative_path,
            status="right_only",
            left_sha256=None,
            right_sha256=right_sha,
            left_size_bytes=None,
            right_size_bytes=right_size,
            detail="artifact appears only in the right package version",
        )
    if right_row is None:
        return PublicationPackageComparisonArtifactRow(
            section=section,
            kind=kind,
            relative_path=relative_path,
            status="left_only",
            left_sha256=left_sha,
            right_sha256=None,
            left_size_bytes=left_size,
            right_size_bytes=None,
            detail="artifact appears only in the left package version",
        )
    if left_sha == right_sha and left_size == right_size:
        return PublicationPackageComparisonArtifactRow(
            section=section,
            kind=kind,
            relative_path=relative_path,
            status="same",
            left_sha256=left_sha,
            right_sha256=right_sha,
            left_size_bytes=left_size,
            right_size_bytes=right_size,
            detail="artifact matches across both package versions",
        )
    detail_parts: list[str] = []
    if left_sha != right_sha:
        detail_parts.append("checksum changed")
    if left_size != right_size:
        detail_parts.append(
            f"size changed from {left_size if left_size is not None else 'missing'} to {right_size if right_size is not None else 'missing'} bytes"
        )
    return PublicationPackageComparisonArtifactRow(
        section=section,
        kind=kind,
        relative_path=relative_path,
        status="changed",
        left_sha256=left_sha,
        right_sha256=right_sha,
        left_size_bytes=left_size,
        right_size_bytes=right_size,
        detail="; ".join(detail_parts),
    )


def load_accession_ids(path: Path) -> set[str]:
    """Load governed accession identifiers from one dataset table."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {
            row["accession"].strip()
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("accession", "").strip()
        }


def load_sequence_ids(path: Path) -> set[str]:
    """Load sequence identifiers from one FASTA dataset surface."""
    return {record.identifier for record in load_fasta_records(path)}


def load_scientific_findings(path: Path) -> dict[str, dict[str, str]]:
    """Load scientific-finding rows keyed by finding id."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {
            row["finding_id"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("finding_id", "").strip()
        }


def inventory_rows_from_manifest(
    package_root: Path,
    manifest: dict[str, object],
) -> list[dict[str, str]]:
    """Resolve and load the artifact inventory referenced by one manifest."""
    inventory_entry = mapping(mapping(manifest, "package_files"), "artifact_inventory")
    inventory_path = package_root / entry_path(inventory_entry)
    return read_tsv_rows(inventory_path)


def package_artifact_rows(
    *,
    left_inventory_rows: list[dict[str, str]],
    right_inventory_rows: list[dict[str, str]],
    left_manifest_path: Path,
    right_manifest_path: Path,
) -> list[PublicationPackageComparisonArtifactRow]:
    """Build artifact comparison rows across both package inventories."""
    left_index = inventory_index(left_inventory_rows)
    right_index = inventory_index(right_inventory_rows)
    relative_paths = sorted(set(left_index) | set(right_index))
    rows = [
        artifact_row(
            relative_path, left_index.get(relative_path), right_index.get(relative_path)
        )
        for relative_path in relative_paths
    ]
    rows.append(
        PublicationPackageComparisonArtifactRow(
            section="package",
            kind="manifest",
            relative_path=left_manifest_path.name,
            status="same"
            if checksum(left_manifest_path) == checksum(right_manifest_path)
            else "changed",
            left_sha256=checksum(left_manifest_path),
            right_sha256=checksum(right_manifest_path),
            left_size_bytes=left_manifest_path.stat().st_size,
            right_size_bytes=right_manifest_path.stat().st_size,
            detail="package manifest checksum matches"
            if checksum(left_manifest_path) == checksum(right_manifest_path)
            else "package manifest checksum changed",
        )
    )
    return rows
