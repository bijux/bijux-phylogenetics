from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class AlignmentRecord:
    """Single FASTA alignment record."""

    identifier: str
    sequence: str


@dataclass(slots=True)
class AlignmentSummary:
    """Summary of an alignment input."""

    path: Path
    sequence_count: int
    alignment_length: int
    min_sequence_length: int
    max_sequence_length: int
    ids: list[str]
    missing_data_fraction: float
    gap_fraction: float
    variable_site_count: int
    parsimony_informative_site_count: int
