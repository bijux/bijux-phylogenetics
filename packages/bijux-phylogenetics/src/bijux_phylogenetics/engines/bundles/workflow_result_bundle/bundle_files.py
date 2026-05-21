from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil

from .contracts import WorkflowResultBundleFile


def copy_bundle_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def record_bundle_file(
    *,
    role: str,
    label: str,
    bundle_root: Path,
    path: Path,
    source_path: Path | None,
) -> WorkflowResultBundleFile:
    relative_path = path.relative_to(bundle_root)
    return WorkflowResultBundleFile(
        role=role,
        label=label,
        relative_path=relative_path,
        sha256=sha256_file(path),
        size_bytes=path.stat().st_size,
        source_path=None if source_path is None else str(source_path),
    )


def maybe_path(value: object) -> Path | None:
    if value is None:
        return None
    return Path(str(value))


def input_label(*, index: int, path: Path) -> str:
    return f"{index:02d}-{path.name}"


def prepared_input_label(path: Path) -> str:
    return f"prepared-{path.name}"


def output_filename(*, label: str, source_path: Path) -> str:
    return (
        f"{label}{source_path.suffix}"
        if source_path.suffix
        else f"{label}-{source_path.name}"
    )


def write_bundle_json(path: Path, payload: dict[str, object]) -> Path:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()
