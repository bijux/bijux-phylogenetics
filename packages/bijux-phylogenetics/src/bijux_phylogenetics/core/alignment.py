from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class AlignmentSummary:
    """Summary of an alignment input."""

    path: Path
    sequence_count: int
    min_sequence_length: int
    max_sequence_length: int

