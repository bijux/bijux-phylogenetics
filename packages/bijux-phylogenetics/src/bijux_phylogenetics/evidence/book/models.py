from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class EvidenceBookValidationIssue:
    path: Path
    message: str


@dataclass(slots=True)
class EvidenceBookValidationReport:
    root: Path
    valid: bool
    issues: list[EvidenceBookValidationIssue]
    bundle_paths: list[Path]
