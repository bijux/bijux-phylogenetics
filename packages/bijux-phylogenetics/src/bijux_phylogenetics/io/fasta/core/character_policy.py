from __future__ import annotations

from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet, AlignmentRecord

_GAP_CHARACTERS = {"-"}
_EXPLICIT_MISSING_CHARACTERS = {"?"}
_DNA_CHARACTERS = set("ACGTNRYSWKMBDHVacgtnryswkmbdhv")
_RNA_CHARACTERS = set("ACGUNRYSWKMBDHVacgunryswkmbdhv")
_NUCLEOTIDE_GC_CHARACTERS = {"G", "C", "g", "c"}
_DNA_CANONICAL = ("A", "C", "G", "T")
_RNA_CANONICAL = ("A", "C", "G", "U")
_PROTEIN_CANONICAL = (
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
)
_PROTEIN_CHARACTERS = set("ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz*")
_DNA_CHARACTERS_UPPER = {character.upper() for character in _DNA_CHARACTERS}
_RNA_CHARACTERS_UPPER = {character.upper() for character in _RNA_CHARACTERS}
_PROTEIN_CHARACTERS_UPPER = {character.upper() for character in _PROTEIN_CHARACTERS}
_SUPPORTED_SEQUENCE_CHARACTERS_UPPER = (
    _DNA_CHARACTERS_UPPER | _RNA_CHARACTERS_UPPER | _PROTEIN_CHARACTERS_UPPER
)
_PROTEIN_EXCLUSIVE_CHARACTERS_UPPER = _PROTEIN_CHARACTERS_UPPER - (
    _DNA_CHARACTERS_UPPER | _RNA_CHARACTERS_UPPER
)
_DNA_AMBIGUITY_UPPER = {"N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}
_RNA_AMBIGUITY_UPPER = {"N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}
_PROTEIN_AMBIGUITY_UPPER = {"B", "J", "X", "Z"}


def _validate_fraction_threshold(threshold: float) -> None:
    if not 0.0 <= threshold <= 1.0:
        raise ValueError(
            f"threshold must be between 0 and 1 inclusive, got {threshold}"
        )


def _is_explicit_missing(residue: str) -> bool:
    return residue in _EXPLICIT_MISSING_CHARACTERS


def _ambiguity_characters_for_alphabet(alphabet: AlignmentAlphabet) -> set[str]:
    if alphabet == "dna":
        return {
            character.lower() for character in _DNA_AMBIGUITY_UPPER
        } | _DNA_AMBIGUITY_UPPER
    if alphabet == "rna":
        return {
            character.lower() for character in _RNA_AMBIGUITY_UPPER
        } | _RNA_AMBIGUITY_UPPER
    if alphabet == "protein":
        return {
            character.lower() for character in _PROTEIN_AMBIGUITY_UPPER
        } | _PROTEIN_AMBIGUITY_UPPER
    return set()


def _is_ambiguity_character(residue: str, *, alphabet: AlignmentAlphabet) -> bool:
    return residue in _ambiguity_characters_for_alphabet(alphabet)


def _is_missing_like(residue: str, *, alphabet: AlignmentAlphabet) -> bool:
    return _is_explicit_missing(residue) or _is_ambiguity_character(
        residue,
        alphabet=alphabet,
    )


def _observed_residues(records: list[AlignmentRecord]) -> list[str]:
    return [
        residue
        for record in records
        for residue in record.sequence
        if residue not in _GAP_CHARACTERS and not _is_explicit_missing(residue)
    ]


def _ordered_sequence_types(types: set[AlignmentAlphabet]) -> list[AlignmentAlphabet]:
    return [alphabet for alphabet in ("dna", "rna", "protein") if alphabet in types]


def _observed_raw_sequence_characters(record: AlignmentRecord) -> set[str]:
    return {
        residue.upper()
        for residue in record.sequence
        if residue not in _GAP_CHARACTERS and not _is_explicit_missing(residue)
    }


def _compatible_raw_sequence_types(record: AlignmentRecord) -> set[AlignmentAlphabet]:
    observed = _observed_raw_sequence_characters(record)
    if not observed:
        return {"dna", "rna", "protein"}
    if observed - _SUPPORTED_SEQUENCE_CHARACTERS_UPPER:
        return set()

    compatible: set[AlignmentAlphabet] = set()
    if observed <= _DNA_CHARACTERS_UPPER:
        compatible.add("dna")
    if observed <= _RNA_CHARACTERS_UPPER:
        compatible.add("rna")
    if observed <= _PROTEIN_CHARACTERS_UPPER:
        compatible.add("protein")
    return compatible


def _allowed_characters_for_sequence_type(
    alphabet: AlignmentAlphabet | None,
) -> set[str] | None:
    if alphabet == "dna":
        return _DNA_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
    if alphabet == "rna":
        return _RNA_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
    if alphabet == "protein":
        return _PROTEIN_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
    return None
