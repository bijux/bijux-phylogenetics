from __future__ import annotations

import hashlib
import json
import shutil
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.errors import EvidenceContractError


@dataclass(slots=True)
class EvidenceFile:
    relative_path: Path
    sha256: str
    size_bytes: int


@dataclass(slots=True)
class EvidenceBundleReport:
    run_root: Path
    output_root: Path
    file_count: int
    files: list[EvidenceFile]


@dataclass(slots=True)
class EvidenceMismatch:
    relative_path: Path
    reason: str


@dataclass(slots=True)
class EvidenceValidationReport:
    bundle_root: Path
    file_count: int
    valid: bool
    mismatches: list[EvidenceMismatch]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bundle_directory(run_root: Path, output_root: Path) -> EvidenceBundleReport:
    """Copy a run directory into a checksummed evidence bundle."""
    if not run_root.exists():
        raise EvidenceContractError(f"run directory not found: {run_root}")
    if not run_root.is_dir():
        raise EvidenceContractError(f"run root is not a directory: {run_root}")

    files_dir = output_root / "files"
    if output_root.exists():
        shutil.rmtree(output_root)
    files_dir.mkdir(parents=True, exist_ok=True)

    bundled_files: list[EvidenceFile] = []
    for source in sorted(path for path in run_root.rglob("*") if path.is_file()):
        relative_path = source.relative_to(run_root)
        destination = files_dir / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        bundled_files.append(
            EvidenceFile(
                relative_path=relative_path,
                sha256=_sha256(source),
                size_bytes=source.stat().st_size,
            )
        )

    manifest_lines = [
        "{",
        f'  "run_root": "{run_root}",',
        '  "files": [',
    ]
    for index, file in enumerate(bundled_files):
        comma = "," if index < len(bundled_files) - 1 else ""
        manifest_lines.extend(
            [
                "    {",
                f'      "relative_path": "{file.relative_path.as_posix()}",',
                f'      "sha256": "{file.sha256}",',
                f'      "size_bytes": {file.size_bytes}',
                f"    }}{comma}",
            ]
        )
    manifest_lines.extend(["  ]", "}"])
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "manifest.json").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    return EvidenceBundleReport(
        run_root=run_root,
        output_root=output_root,
        file_count=len(bundled_files),
        files=bundled_files,
    )


def validate_bundle(bundle_root: Path) -> EvidenceValidationReport:
    """Validate an evidence bundle against its manifest checksums."""
    manifest_path = bundle_root / "manifest.json"
    files_root = bundle_root / "files"
    if not manifest_path.exists():
        raise EvidenceContractError(f"evidence manifest not found: {manifest_path}")
    if not files_root.exists():
        raise EvidenceContractError(f"evidence files directory not found: {files_root}")

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_files = payload.get("files")
    if not isinstance(manifest_files, list):
        raise EvidenceContractError(f"evidence manifest is missing a files list: {manifest_path}")

    mismatches: list[EvidenceMismatch] = []
    for entry in manifest_files:
        relative_path = Path(str(entry["relative_path"]))
        expected_sha = str(entry["sha256"])
        expected_size = int(entry["size_bytes"])
        candidate = files_root / relative_path
        if not candidate.exists():
            mismatches.append(EvidenceMismatch(relative_path=relative_path, reason="missing_file"))
            continue
        if candidate.stat().st_size != expected_size:
            mismatches.append(EvidenceMismatch(relative_path=relative_path, reason="size_mismatch"))
            continue
        if _sha256(candidate) != expected_sha:
            mismatches.append(EvidenceMismatch(relative_path=relative_path, reason="checksum_mismatch"))

    return EvidenceValidationReport(
        bundle_root=bundle_root,
        file_count=len(manifest_files),
        valid=not mismatches,
        mismatches=mismatches,
    )
