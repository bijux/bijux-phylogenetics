from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class EvidenceEntry:
    """Single file tracked in an evidence bundle."""

    relative_path: Path
    sha256: str
    size_bytes: int

