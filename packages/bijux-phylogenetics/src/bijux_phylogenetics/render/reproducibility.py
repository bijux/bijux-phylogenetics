from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class FigureReproducibilityArtifact:
    """One reproducibility-tracked figure input or output artifact."""

    label: str
    path: str
    checksum: str


@dataclass(frozen=True, slots=True)
class FigureReproducibilityFilter:
    """One explicit filtering or selection rule applied before figure rendering."""

    name: str
    value: str
    detail: str


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_ready(payload: object) -> object:
    return json.loads(json.dumps(payload, default=str))


def _artifact_rows(
    artifacts: list[tuple[str, Path]] | None,
) -> list[FigureReproducibilityArtifact]:
    if artifacts is None:
        return []
    return [
        FigureReproducibilityArtifact(
            label=label,
            path=str(path),
            checksum=_checksum(path),
        )
        for label, path in artifacts
    ]


def build_figure_reproducibility_manifest(
    *,
    report_kind: str,
    input_files: list[tuple[str, Path]],
    generated_figures: list[tuple[str, Path]],
    generated_tables: list[tuple[str, Path]],
    filters: list[FigureReproducibilityFilter] | None,
    model: dict[str, Any],
    settings: dict[str, Any],
    linked_artifacts: list[tuple[str, Path]] | None = None,
) -> dict[str, object]:
    """Build one standardized reproducibility manifest for a figure package."""
    return {
        "report_kind": report_kind,
        "format_version": 1,
        "input_files": [asdict(row) for row in _artifact_rows(input_files)],
        "filters": [] if filters is None else [asdict(row) for row in filters],
        "model": _json_ready(model),
        "settings": _json_ready(settings),
        "generated_figures": [asdict(row) for row in _artifact_rows(generated_figures)],
        "generated_tables": [asdict(row) for row in _artifact_rows(generated_tables)],
        "linked_artifacts": [asdict(row) for row in _artifact_rows(linked_artifacts)],
    }


def write_figure_reproducibility_manifest(
    path: Path,
    *,
    report_kind: str,
    input_files: list[tuple[str, Path]],
    generated_figures: list[tuple[str, Path]],
    generated_tables: list[tuple[str, Path]],
    filters: list[FigureReproducibilityFilter] | None,
    model: dict[str, Any],
    settings: dict[str, Any],
    linked_artifacts: list[tuple[str, Path]] | None = None,
) -> dict[str, object]:
    """Write one standardized reproducibility manifest for a figure package."""
    manifest = build_figure_reproducibility_manifest(
        report_kind=report_kind,
        input_files=input_files,
        generated_figures=generated_figures,
        generated_tables=generated_tables,
        filters=filters,
        model=model,
        settings=settings,
        linked_artifacts=linked_artifacts,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest
