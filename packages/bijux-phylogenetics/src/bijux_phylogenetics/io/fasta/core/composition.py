from __future__ import annotations

from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    SequenceGCContent,
)

from .character_policy import (
    _DNA_CANONICAL,
    _DNA_CHARACTERS,
    _NUCLEOTIDE_GC_CHARACTERS,
    _PROTEIN_CANONICAL,
    _PROTEIN_CHARACTERS,
    _RNA_CANONICAL,
    _RNA_CHARACTERS,
    _observed_residues,
)


def infer_alignment_alphabet(records: list[AlignmentRecord]) -> AlignmentAlphabet:
    """Infer whether an alignment is DNA, RNA, protein, or unknown."""
    observed = _observed_residues(records)
    if not observed:
        return "unknown"
    characters = set(observed)
    if characters <= _DNA_CHARACTERS and "U" not in {
        residue.upper() for residue in observed
    }:
        return "dna"
    if characters <= _RNA_CHARACTERS and "T" not in {
        residue.upper() for residue in observed
    }:
        return "rna"
    if characters <= _PROTEIN_CHARACTERS:
        return "protein"
    return "unknown"


def _normalized_frequency(
    values: list[str],
    alphabet: tuple[str, ...],
) -> dict[str, float]:
    if not values:
        return {}
    total = len(values)
    return {
        character: round(
            sum(1 for value in values if value.upper() == character) / total,
            15,
        )
        for character in alphabet
        if any(value.upper() == character for value in values)
    }


def compute_nucleotide_composition(
    records: list[AlignmentRecord],
    *,
    alphabet: AlignmentAlphabet,
) -> dict[str, float]:
    """Compute canonical nucleotide composition for DNA or RNA alignments."""
    observed = _observed_residues(records)
    if alphabet == "dna":
        return _normalized_frequency(observed, _DNA_CANONICAL)
    if alphabet == "rna":
        return _normalized_frequency(observed, _RNA_CANONICAL)
    return {}


def compute_amino_acid_composition(
    records: list[AlignmentRecord],
    *,
    alphabet: AlignmentAlphabet,
) -> dict[str, float]:
    """Compute canonical amino-acid composition for protein alignments."""
    if alphabet != "protein":
        return {}
    return _normalized_frequency(_observed_residues(records), _PROTEIN_CANONICAL)


def _sequence_gc_fraction(sequence: str) -> float | None:
    comparable = [
        residue for residue in sequence if residue.upper() in {"A", "C", "G", "T", "U"}
    ]
    if not comparable:
        return None
    return round(
        sum(1 for residue in comparable if residue in _NUCLEOTIDE_GC_CHARACTERS)
        / len(comparable),
        15,
    )


def compute_per_sequence_gc_content(
    records: list[AlignmentRecord],
    *,
    alphabet: AlignmentAlphabet,
) -> list[SequenceGCContent]:
    """Compute GC content for each sequence when the alignment is nucleotide-like."""
    if alphabet not in {"dna", "rna"}:
        return []
    return [
        SequenceGCContent(
            identifier=record.identifier,
            gc_fraction=_sequence_gc_fraction(record.sequence),
        )
        for record in records
    ]


def compute_whole_alignment_gc_content(
    records: list[AlignmentRecord],
    *,
    alphabet: AlignmentAlphabet,
) -> float | None:
    """Compute whole-alignment GC content for nucleotide alignments."""
    if alphabet not in {"dna", "rna"}:
        return None
    comparable = [
        residue
        for record in records
        for residue in record.sequence
        if residue.upper() in {"A", "C", "G", "T", "U"}
    ]
    if not comparable:
        return None
    return round(
        sum(1 for residue in comparable if residue in _NUCLEOTIDE_GC_CHARACTERS)
        / len(comparable),
        15,
    )
