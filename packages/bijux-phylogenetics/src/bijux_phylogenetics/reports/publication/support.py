from __future__ import annotations

import csv
from hashlib import sha256
import json
from pathlib import Path

SUPPORTED_PUBLICATION_PACKAGE_KIND = "rabies_cross_host_geography_package"


def checksum(path: Path) -> str:
    """Return one sha256 checksum for a stored publication artifact."""

    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def artifact_kind(relative_path: str) -> str:
    """Classify one package artifact by durable output kind."""

    path = Path(relative_path)
    suffix = path.suffix.lower()
    if path.name.endswith(".manifest.json"):
        return "manifest"
    if suffix in {".html", ".htm"}:
        return "report"
    if suffix == ".md":
        return "markdown"
    if suffix == ".json":
        return "json"
    if suffix == ".tsv":
        return "table"
    if suffix == ".svg":
        return "figure"
    if suffix == ".log":
        return "log"
    if suffix in {".nwk", ".tree"}:
        return "tree"
    if suffix in {".aln", ".fasta"}:
        return "alignment"
    if suffix == ".csv":
        return "metadata"
    return "artifact"


def mapping(payload: dict[str, object], key: str) -> dict[str, object]:
    """Return one mapping value or an empty mapping when the key is absent."""

    value = payload.get(key)
    if isinstance(value, dict):
        return value
    return {}


def text(payload: object) -> str:
    """Normalize one manifest value into trimmed text."""

    return str(payload).strip()


def entry_path(payload: dict[str, object]) -> str:
    """Return one stored relative path from a manifest entry."""

    value = payload.get("path")
    return value.strip() if isinstance(value, str) else ""


def entry_checksum(payload: dict[str, object]) -> str | None:
    """Return one stored checksum from a manifest entry when present."""

    value = payload.get("checksum")
    if isinstance(value, str) and value.strip():
        return value.strip()
    return None


def entry_size(payload: dict[str, object]) -> int | None:
    """Return one stored size from a manifest or inventory entry when present."""

    value = payload.get("size_bytes")
    if value is None:
        return None
    normalized = str(value).strip()
    if not normalized:
        return None
    return int(normalized)


def read_manifest(path: Path) -> dict[str, object]:
    """Read one stored publication package manifest."""

    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} does not contain one JSON object")
    return payload


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    """Read one package TSV artifact as tab-delimited rows."""

    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def manifest_file_entries(
    manifest: dict[str, object],
) -> list[tuple[str, str, str, str | None]]:
    """Return manifest-declared package, dataset, and workflow files."""

    entries: list[tuple[str, str, str, str | None]] = []
    for block_name in ("package_files", "dataset_files", "workflow_files"):
        block = mapping(manifest, block_name)
        for entry_name, entry_payload in block.items():
            if not isinstance(entry_payload, dict):
                continue
            relative_path = entry_path(entry_payload)
            if not relative_path:
                continue
            entries.append(
                (
                    block_name,
                    entry_name,
                    relative_path,
                    entry_checksum(entry_payload),
                )
            )
    return entries


def section_counts(rows: list[dict[str, str]]) -> dict[str, int]:
    """Count inventory rows by declared section."""

    counts: dict[str, int] = {}
    for row in rows:
        section = row.get("section", "").strip()
        if not section:
            continue
        counts[section] = counts.get(section, 0) + 1
    return counts


def ignored_package_prefixes(report_kind: str) -> tuple[str, ...]:
    """Return package-relative prefixes that are governed reference content."""

    if report_kind == SUPPORTED_PUBLICATION_PACKAGE_KIND:
        return ("dataset/expected/",)
    return ()
