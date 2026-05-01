from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import UTC, datetime
import hashlib
from importlib import metadata
import json
from pathlib import Path
import platform
import sys


@dataclass(slots=True)
class RunManifest:
    """Reproducibility manifest for one command execution."""

    command: str
    arguments: list[str]
    input_paths: list[str]
    output_paths: list[str]
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    package_version: str
    python_version: str
    dependency_versions: dict[str, str]
    timestamp_utc: str
    host_platform: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _dependency_version(name: str) -> str:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return "unavailable"


def build_run_manifest(
    *,
    command: str,
    arguments: list[str],
    input_paths: list[Path | str],
    output_paths: list[Path | str] | None = None,
) -> RunManifest:
    """Build a deterministic manifest for one command execution."""
    normalized_input_paths = [Path(item) for item in input_paths]
    normalized_output_paths = [Path(item) for item in (output_paths or [])]
    return RunManifest(
        command=command,
        arguments=list(arguments),
        input_paths=[str(path) for path in normalized_input_paths],
        output_paths=[str(path) for path in normalized_output_paths],
        input_checksums={
            str(path): _sha256(path)
            for path in normalized_input_paths
            if path.exists() and path.is_file()
        },
        output_checksums={
            str(path): _sha256(path)
            for path in normalized_output_paths
            if path.exists() and path.is_file()
        },
        package_version=_dependency_version("bijux-phylogenetics"),
        python_version=sys.version.split()[0],
        dependency_versions={
            "biopython": _dependency_version("biopython"),
            "bijux-phylogenetics": _dependency_version("bijux-phylogenetics"),
        },
        timestamp_utc=datetime.now(UTC).isoformat(),
        host_platform=platform.platform(),
    )


def write_run_manifest(path: Path, manifest: RunManifest) -> Path:
    """Serialize a run manifest to JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(manifest), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path
