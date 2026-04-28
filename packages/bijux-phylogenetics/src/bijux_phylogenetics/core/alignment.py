from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class AlignmentRecord:
    """Single FASTA alignment record."""

    identifier: str
    sequence: str


@dataclass(frozen=True, slots=True)
class SequenceMissingness:
    """Missing-data fraction for one alignment sequence."""

    identifier: str
    missing_fraction: float


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
    per_sequence_missingness: list[SequenceMissingness]
    constant_site_count: int
    variable_site_count: int
    parsimony_informative_site_count: int


@dataclass(slots=True)
class AlignmentLinkageReport:
    """Summary of how an alignment links against a tree tip set."""

    tree_path: Path
    alignment_path: Path
    tree_taxa: int
    alignment_ids: int
    linked_taxa: int
    usable_taxa: list[str]
    missing_from_alignment: list[str]
    extra_alignment_ids: list[str]
