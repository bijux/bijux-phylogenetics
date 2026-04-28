from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


AlignmentAlphabet = str


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


@dataclass(frozen=True, slots=True)
class SiteMissingness:
    """Missing-data fraction for one alignment column."""

    position: int
    missing_fraction: float


@dataclass(frozen=True, slots=True)
class InvalidAlignmentCharacter:
    """One sequence character invalid for a declared alignment alphabet."""

    identifier: str
    position: int
    character: str


@dataclass(frozen=True, slots=True)
class SequenceGCContent:
    """GC content summary for one sequence."""

    identifier: str
    gc_fraction: float | None


@dataclass(frozen=True, slots=True)
class SequenceCompositionOutlier:
    """One sequence whose composition deviates strongly from the alignment baseline."""

    identifier: str
    deviation: float


@dataclass(frozen=True, slots=True)
class DuplicateSequenceGroup:
    """Identifiers sharing the exact same aligned sequence."""

    identifiers: list[str]
    sequence: str


@dataclass(frozen=True, slots=True)
class NearDuplicateSequencePair:
    """Pair of sequences above a caller-provided identity threshold."""

    left_identifier: str
    right_identifier: str
    identity: float
    comparable_sites: int


@dataclass(slots=True)
class AlignmentQualityReport:
    """Higher-level alignment quality report built from composition and identity diagnostics."""

    path: Path
    sequence_count: int
    alignment_length: int
    missing_data_fraction: float
    gap_fraction: float
    variable_site_count: int
    parsimony_informative_site_count: int
    inferred_alphabet: AlignmentAlphabet
    invalid_characters: list[InvalidAlignmentCharacter]
    composition_outliers: list[SequenceCompositionOutlier]
    duplicate_sequence_groups: list[DuplicateSequenceGroup]
    near_duplicate_pairs: list[NearDuplicateSequencePair]
    warnings: list[str]


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
    per_site_missingness: list[SiteMissingness]
    all_gap_columns: list[int]
    all_missing_columns: list[int]
    constant_site_count: int
    variable_site_count: int
    parsimony_informative_site_count: int
    inferred_alphabet: AlignmentAlphabet
    invalid_characters: list[InvalidAlignmentCharacter]
    nucleotide_composition: dict[str, float]
    amino_acid_composition: dict[str, float]
    per_sequence_gc_content: list[SequenceGCContent]
    whole_alignment_gc_content: float | None
    composition_outliers: list[SequenceCompositionOutlier]
    duplicate_sequence_groups: list[DuplicateSequenceGroup]
    near_duplicate_pairs: list[NearDuplicateSequencePair]


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
