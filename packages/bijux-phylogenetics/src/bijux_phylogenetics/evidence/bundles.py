from __future__ import annotations

import hashlib
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
