from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class ManifestReplayDrift:
    kind: str
    label: str
    expected: str
    observed: str | None
    matched: bool


@dataclass(slots=True)
class ManifestReplayComparison:
    label: str
    status: str
    detail: str


@dataclass(slots=True)
class ManifestReplayReport:
    manifest_path: Path
    workflow: str
    replay_out_dir: Path
    replay_manifest_path: Path
    input_drift: list[ManifestReplayDrift]
    engine_version_drift: list[ManifestReplayDrift]
    comparisons: list[ManifestReplayComparison]
    input_drift_detected: bool
    engine_version_drift_detected: bool
    outputs_equivalent: bool
    notes: list[str]
