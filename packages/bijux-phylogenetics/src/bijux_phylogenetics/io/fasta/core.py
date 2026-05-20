from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from statistics import median

from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    FastaDuplicateIdentifier,
    FastaEmptySequence,
    FastaIllegalCharacter,
    FastaInputSummary,
    FastaSequenceTypeReport,
    InvalidAlignmentCharacter,
    SequenceGCContent,
    SequenceLengthOutlier,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

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


@dataclass(slots=True)
class _RawFastaScan:
    path: Path
    record_count: int
    lengths_by_identifier: list[tuple[str, int]]
    total_residue_count: int
    duplicate_positions: dict[str, list[int]]
    empty_sequences: list[FastaEmptySequence]
    illegal_supported: list[FastaIllegalCharacter]
    illegal_dna: list[FastaIllegalCharacter]
    illegal_rna: list[FastaIllegalCharacter]
    illegal_protein: list[FastaIllegalCharacter]
    thymine_record_count: int
    uracil_record_count: int
    protein_signal_record_count: int
    invalid_record_count: int
    shared_compatible: set[AlignmentAlphabet] | None


def _validate_fraction_threshold(threshold: float) -> None:
    if not 0.0 <= threshold <= 1.0:
        raise ValueError(
            f"threshold must be between 0 and 1 inclusive, got {threshold}"
        )


def _median_absolute_deviation(values: list[float]) -> float:
    if not values:
        return 0.0
    center = median(values)
    return float(median([abs(value - center) for value in values]))


def _robust_z_score(value: float, values: list[float]) -> float | None:
    mad = _median_absolute_deviation(values)
    if mad == 0.0:
        return None
    return round(0.6744897501960817 * (value - float(median(values))) / mad, 15)


def _is_explicit_missing(residue: str) -> bool:
    return residue in _EXPLICIT_MISSING_CHARACTERS


def _ambiguity_characters_for_alphabet(alphabet: AlignmentAlphabet) -> set[str]:
    if alphabet == "dna":
        return {character.lower() for character in _DNA_AMBIGUITY_UPPER} | _DNA_AMBIGUITY_UPPER
    if alphabet == "rna":
        return {character.lower() for character in _RNA_AMBIGUITY_UPPER} | _RNA_AMBIGUITY_UPPER
    if alphabet == "protein":
        return {character.lower() for character in _PROTEIN_AMBIGUITY_UPPER} | _PROTEIN_AMBIGUITY_UPPER
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


def _detect_duplicate_identifiers(
    records: list[AlignmentRecord],
) -> list[FastaDuplicateIdentifier]:
    positions: dict[str, list[int]] = {}
    for index, record in enumerate(records, start=1):
        positions.setdefault(record.identifier, []).append(index)
    return [
        FastaDuplicateIdentifier(
            identifier=identifier,
            occurrences=len(record_indices),
            record_indices=record_indices,
        )
        for identifier, record_indices in sorted(positions.items())
        if len(record_indices) > 1
    ]


def _detect_empty_sequences(records: list[AlignmentRecord]) -> list[FastaEmptySequence]:
    return [
        FastaEmptySequence(identifier=record.identifier, record_index=index)
        for index, record in enumerate(records, start=1)
        if not record.sequence
    ]


def _detect_illegal_sequence_characters(
    records: list[AlignmentRecord],
    *,
    alphabet: AlignmentAlphabet | None,
) -> list[FastaIllegalCharacter]:
    allowed = _allowed_characters_for_sequence_type(alphabet)
    if allowed is None:
        allowed = (
            _DNA_CHARACTERS
            | _RNA_CHARACTERS
            | _PROTEIN_CHARACTERS
            | _GAP_CHARACTERS
            | _EXPLICIT_MISSING_CHARACTERS
        )
    return [
        FastaIllegalCharacter(
            identifier=record.identifier,
            record_index=record_index,
            position=position,
            character=residue,
        )
        for record_index, record in enumerate(records, start=1)
        for position, residue in enumerate(record.sequence, start=1)
        if residue not in allowed
    ]


def _normalize_fasta_identifier(identifier: str) -> str:
    normalized = [
        character if character.isalnum() or character in {".", "_", "-"} else "_"
        for character in identifier.strip()
    ]
    collapsed = "".join(normalized)
    while "__" in collapsed:
        collapsed = collapsed.replace("__", "_")
    collapsed = collapsed.strip("._-")
    return collapsed or "sequence"


def _process_raw_fasta_record(
    *,
    identifier: str,
    sequence: str,
    record_index: int,
    lengths_by_identifier: list[tuple[str, int]],
    duplicate_positions: dict[str, list[int]],
    empty_sequences: list[FastaEmptySequence],
    illegal_supported: list[FastaIllegalCharacter],
    illegal_dna: list[FastaIllegalCharacter],
    illegal_rna: list[FastaIllegalCharacter],
    illegal_protein: list[FastaIllegalCharacter],
    thymine_record_count_ref: list[int],
    uracil_record_count_ref: list[int],
    protein_signal_record_count_ref: list[int],
    invalid_record_count_ref: list[int],
    shared_compatible_ref: list[set[AlignmentAlphabet] | None],
) -> int:
    record = AlignmentRecord(identifier=identifier, sequence=sequence)
    lengths_by_identifier.append((identifier, len(sequence)))
    duplicate_positions.setdefault(identifier, []).append(record_index)
    if not sequence:
        empty_sequences.append(
            FastaEmptySequence(identifier=identifier, record_index=record_index)
        )
    observed = _observed_raw_sequence_characters(record)
    compatible = _compatible_raw_sequence_types(record)
    if compatible:
        if shared_compatible_ref[0] is None:
            shared_compatible_ref[0] = set(compatible)
        else:
            shared_compatible_ref[0] &= compatible
    else:
        invalid_record_count_ref[0] += 1
    if "T" in observed:
        thymine_record_count_ref[0] += 1
    if "U" in observed:
        uracil_record_count_ref[0] += 1
    if observed & _PROTEIN_EXCLUSIVE_CHARACTERS_UPPER:
        protein_signal_record_count_ref[0] += 1
    for position, residue in enumerate(sequence, start=1):
        row = FastaIllegalCharacter(
            identifier=identifier,
            record_index=record_index,
            position=position,
            character=residue,
        )
        if residue not in (
            _DNA_CHARACTERS
            | _RNA_CHARACTERS
            | _PROTEIN_CHARACTERS
            | _GAP_CHARACTERS
            | _EXPLICIT_MISSING_CHARACTERS
        ):
            illegal_supported.append(row)
        if residue not in (
            _DNA_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
        ):
            illegal_dna.append(row)
        if residue not in (
            _RNA_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
        ):
            illegal_rna.append(row)
        if residue not in (
            _PROTEIN_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
        ):
            illegal_protein.append(row)
    return len(sequence)


def _scan_raw_fasta(path: Path) -> _RawFastaScan:
    if not path.exists():
        raise FileNotFoundError(f"alignment file not found: {path}")

    record_count = 0
    lengths_by_identifier: list[tuple[str, int]] = []
    total_residue_count = 0
    duplicate_positions: dict[str, list[int]] = {}
    empty_sequences: list[FastaEmptySequence] = []
    illegal_supported: list[FastaIllegalCharacter] = []
    illegal_dna: list[FastaIllegalCharacter] = []
    illegal_rna: list[FastaIllegalCharacter] = []
    illegal_protein: list[FastaIllegalCharacter] = []
    thymine_record_count_ref = [0]
    uracil_record_count_ref = [0]
    protein_signal_record_count_ref = [0]
    invalid_record_count_ref = [0]
    shared_compatible_ref: list[set[AlignmentAlphabet] | None] = [None]
    current_identifier: str | None = None
    current_sequence: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_identifier is not None:
                    record_count += 1
                    total_residue_count += _process_raw_fasta_record(
                        identifier=current_identifier,
                        sequence="".join(current_sequence),
                        record_index=record_count,
                        lengths_by_identifier=lengths_by_identifier,
                        duplicate_positions=duplicate_positions,
                        empty_sequences=empty_sequences,
                        illegal_supported=illegal_supported,
                        illegal_dna=illegal_dna,
                        illegal_rna=illegal_rna,
                        illegal_protein=illegal_protein,
                        thymine_record_count_ref=thymine_record_count_ref,
                        uracil_record_count_ref=uracil_record_count_ref,
                        protein_signal_record_count_ref=protein_signal_record_count_ref,
                        invalid_record_count_ref=invalid_record_count_ref,
                        shared_compatible_ref=shared_compatible_ref,
                    )
                current_identifier = line[1:].strip()
                current_sequence = []
                continue
            if current_identifier is None:
                raise InvalidAlignmentError(
                    f"alignment sequence appears before any FASTA header in {path}"
                )
            current_sequence.append(line)

    if current_identifier is not None:
        record_count += 1
        total_residue_count += _process_raw_fasta_record(
            identifier=current_identifier,
            sequence="".join(current_sequence),
            record_index=record_count,
            lengths_by_identifier=lengths_by_identifier,
            duplicate_positions=duplicate_positions,
            empty_sequences=empty_sequences,
            illegal_supported=illegal_supported,
            illegal_dna=illegal_dna,
            illegal_rna=illegal_rna,
            illegal_protein=illegal_protein,
            thymine_record_count_ref=thymine_record_count_ref,
            uracil_record_count_ref=uracil_record_count_ref,
            protein_signal_record_count_ref=protein_signal_record_count_ref,
            invalid_record_count_ref=invalid_record_count_ref,
            shared_compatible_ref=shared_compatible_ref,
        )

    if record_count == 0:
        raise InvalidAlignmentError(f"alignment contains no FASTA records: {path}")

    return _RawFastaScan(
        path=path,
        record_count=record_count,
        lengths_by_identifier=lengths_by_identifier,
        total_residue_count=total_residue_count,
        duplicate_positions=duplicate_positions,
        empty_sequences=empty_sequences,
        illegal_supported=illegal_supported,
        illegal_dna=illegal_dna,
        illegal_rna=illegal_rna,
        illegal_protein=illegal_protein,
        thymine_record_count=thymine_record_count_ref[0],
        uracil_record_count=uracil_record_count_ref[0],
        protein_signal_record_count=protein_signal_record_count_ref[0],
        invalid_record_count=invalid_record_count_ref[0],
        shared_compatible=shared_compatible_ref[0],
    )


def _build_fasta_sequence_type_report(
    *,
    path: Path,
    record_count: int,
    shared_compatible: set[AlignmentAlphabet] | None,
    thymine_record_count: int,
    uracil_record_count: int,
    protein_signal_record_count: int,
    invalid_record_count: int,
) -> FastaSequenceTypeReport:
    warnings: list[str] = []
    compatible_types = (
        [] if shared_compatible is None else _ordered_sequence_types(shared_compatible)
    )
    if invalid_record_count:
        warnings.append("input contains unsupported sequence characters")
        return FastaSequenceTypeReport(
            path=path,
            record_count=record_count,
            detected_type="invalid",
            selected_type=None,
            compatible_types=[],
            confidence="blocked",
            thymine_record_count=thymine_record_count,
            uracil_record_count=uracil_record_count,
            protein_signal_record_count=protein_signal_record_count,
            invalid_record_count=invalid_record_count,
            note=(
                "one or more records contain characters outside the supported DNA, "
                "RNA, and protein alphabets"
            ),
            warnings=warnings,
        )

    if not compatible_types:
        warnings.append("input contains conflicting sequence-type signals")
        if thymine_record_count and uracil_record_count:
            note = (
                "records mix thymine-bearing and uracil-bearing sequences, so one "
                "shared nucleotide workflow is not defensible"
            )
        else:
            note = (
                "records do not share one compatible DNA, RNA, or protein type "
                "without an explicit caller override"
            )
        return FastaSequenceTypeReport(
            path=path,
            record_count=record_count,
            detected_type="mixed",
            selected_type=None,
            compatible_types=[],
            confidence="blocked",
            thymine_record_count=thymine_record_count,
            uracil_record_count=uracil_record_count,
            protein_signal_record_count=protein_signal_record_count,
            invalid_record_count=0,
            note=note,
            warnings=warnings,
        )

    if protein_signal_record_count:
        return FastaSequenceTypeReport(
            path=path,
            record_count=record_count,
            detected_type="protein",
            selected_type="protein",
            compatible_types=compatible_types,
            confidence="high",
            thymine_record_count=thymine_record_count,
            uracil_record_count=uracil_record_count,
            protein_signal_record_count=protein_signal_record_count,
            invalid_record_count=0,
            note="protein-exclusive residues were observed in the raw FASTA input",
            warnings=warnings,
        )

    if uracil_record_count and not thymine_record_count:
        return FastaSequenceTypeReport(
            path=path,
            record_count=record_count,
            detected_type="rna",
            selected_type="rna",
            compatible_types=compatible_types,
            confidence="high",
            thymine_record_count=thymine_record_count,
            uracil_record_count=uracil_record_count,
            protein_signal_record_count=protein_signal_record_count,
            invalid_record_count=0,
            note="uracil was observed without thymine, which strongly supports RNA",
            warnings=warnings,
        )

    if "dna" in compatible_types:
        confidence = "medium" if thymine_record_count else "low"
        if confidence == "medium":
            warnings.append(
                "automatic sequence type defaults to dna from nucleotide-like characters that remain protein-compatible by alphabet alone"
            )
            note = (
                "thymine was observed and no protein-exclusive residues were found, "
                "so the raw input defaults to DNA"
            )
        else:
            warnings.append(
                "automatic sequence type defaults to dna from characters shared by DNA, RNA, and protein alphabets"
            )
            note = (
                "the raw input uses only characters shared across multiple "
                "alphabets, so DNA is chosen as the default but an explicit "
                "sequence type remains safer when the biological context is unclear"
            )
        return FastaSequenceTypeReport(
            path=path,
            record_count=record_count,
            detected_type="dna",
            selected_type="dna",
            compatible_types=compatible_types,
            confidence=confidence,
            thymine_record_count=thymine_record_count,
            uracil_record_count=uracil_record_count,
            protein_signal_record_count=protein_signal_record_count,
            invalid_record_count=0,
            note=note,
            warnings=warnings,
        )

    return FastaSequenceTypeReport(
        path=path,
        record_count=record_count,
        detected_type="unknown",
        selected_type=None,
        compatible_types=compatible_types,
        confidence="blocked",
        thymine_record_count=thymine_record_count,
        uracil_record_count=uracil_record_count,
        protein_signal_record_count=protein_signal_record_count,
        invalid_record_count=0,
        note="the raw input could not be resolved to a supported sequence type",
        warnings=["input sequence type could not be inferred confidently"],
    )


def _build_fasta_input_summary_from_scan(
    scan: _RawFastaScan,
    *,
    sequence_type: AlignmentAlphabet | None,
    sequence_type_report: FastaSequenceTypeReport,
) -> FastaInputSummary:
    lengths = [length for _, length in scan.lengths_by_identifier]
    inferred_alphabet = (
        sequence_type
        if sequence_type is not None
        else (
            "unknown"
            if sequence_type_report.selected_type is None
            else sequence_type_report.selected_type
        )
    )
    return FastaInputSummary(
        path=scan.path,
        sequence_count=scan.record_count,
        unique_identifier_count=len(scan.duplicate_positions),
        empty_sequence_count=len(scan.empty_sequences),
        min_sequence_length=min(lengths),
        max_sequence_length=max(lengths),
        median_sequence_length=float(median(lengths)),
        total_residue_count=scan.total_residue_count,
        inferred_alphabet=inferred_alphabet,
    )


def detect_fasta_sequence_type(
    path: Path,
    *,
    records: list[AlignmentRecord] | None = None,
) -> FastaSequenceTypeReport:
    """Classify raw FASTA records as DNA, RNA, protein, mixed, or invalid."""
    if records is None:
        scan = _scan_raw_fasta(path)
        return _build_fasta_sequence_type_report(
            path=path,
            record_count=scan.record_count,
            shared_compatible=scan.shared_compatible,
            thymine_record_count=scan.thymine_record_count,
            uracil_record_count=scan.uracil_record_count,
            protein_signal_record_count=scan.protein_signal_record_count,
            invalid_record_count=scan.invalid_record_count,
        )
    shared_compatible: set[AlignmentAlphabet] | None = None
    thymine_record_count = 0
    uracil_record_count = 0
    protein_signal_record_count = 0
    invalid_record_count = 0
    for record in records:
        observed = _observed_raw_sequence_characters(record)
        compatible = _compatible_raw_sequence_types(record)
        if compatible:
            if shared_compatible is None:
                shared_compatible = set(compatible)
            else:
                shared_compatible &= compatible
        else:
            invalid_record_count += 1
        if "T" in observed:
            thymine_record_count += 1
        if "U" in observed:
            uracil_record_count += 1
        if observed & _PROTEIN_EXCLUSIVE_CHARACTERS_UPPER:
            protein_signal_record_count += 1
    return _build_fasta_sequence_type_report(
        path=path,
        record_count=len(records),
        shared_compatible=shared_compatible,
        thymine_record_count=thymine_record_count,
        uracil_record_count=uracil_record_count,
        protein_signal_record_count=protein_signal_record_count,
        invalid_record_count=invalid_record_count,
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


def _detect_invalid_alignment_characters_records(
    records: list[AlignmentRecord],
    *,
    alphabet: AlignmentAlphabet,
) -> list[InvalidAlignmentCharacter]:
    allowed = _allowed_characters_for_sequence_type(alphabet)
    if allowed is None:
        raise ValueError(f"unsupported declared alphabet: {alphabet}")

    invalid: list[InvalidAlignmentCharacter] = []
    for record in records:
        for position, residue in enumerate(record.sequence, start=1):
            if residue not in allowed:
                invalid.append(
                    InvalidAlignmentCharacter(
                        identifier=record.identifier,
                        position=position,
                        character=residue,
                    )
                )
    return invalid


def detect_invalid_alignment_characters(
    path: Path,
    *,
    alphabet: AlignmentAlphabet,
) -> list[InvalidAlignmentCharacter]:
    """List sequence characters invalid for the declared alphabet."""
    records = load_fasta_alignment(path)
    return _detect_invalid_alignment_characters_records(records, alphabet=alphabet)


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
        residue
        for residue in sequence
        if residue.upper() in {"A", "C", "G", "T", "U"}
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


def load_permissive_fasta_records(path: Path) -> list[AlignmentRecord]:
    """Load FASTA records while preserving duplicates and empty sequences."""
    if not path.exists():
        raise FileNotFoundError(f"alignment file not found: {path}")

    records: list[AlignmentRecord] = []
    current_identifier: str | None = None
    current_sequence: list[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_identifier is not None:
                    records.append(
                        AlignmentRecord(
                            identifier=current_identifier,
                            sequence="".join(current_sequence),
                        )
                    )
                current_identifier = line[1:].strip()
                current_sequence = []
                continue
            if current_identifier is None:
                raise InvalidAlignmentError(
                    f"alignment sequence appears before any FASTA header in {path}"
                )
            current_sequence.append(line)

    if current_identifier is not None:
        records.append(
            AlignmentRecord(
                identifier=current_identifier,
                sequence="".join(current_sequence),
            )
        )

    if not records:
        raise InvalidAlignmentError(f"alignment contains no FASTA records: {path}")

    return records


def load_fasta_records(path: Path) -> list[AlignmentRecord]:
    """Load FASTA records without imposing equal-length alignment requirements."""
    records = load_permissive_fasta_records(path)
    duplicate_rows = _detect_duplicate_identifiers(records)
    if duplicate_rows:
        raise InvalidAlignmentError(
            "alignment contains duplicate sequence ids: "
            + ", ".join(row.identifier for row in duplicate_rows)
        )
    return records


def load_fasta_alignment(path: Path) -> list[AlignmentRecord]:
    """Load FASTA records using the repository's deterministic alignment contract."""
    records = load_fasta_records(path)
    lengths = [len(record.sequence) for record in records]
    if min(lengths) != max(lengths):
        raise InvalidAlignmentError(
            f"alignment contains unequal sequence lengths: min={min(lengths)} max={max(lengths)}"
        )
    return records


def write_fasta_alignment(path: Path, records: list[AlignmentRecord]) -> Path:
    """Write FASTA records using a deterministic single-line sequence layout."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for record in records:
        lines.append(f">{record.identifier}")
        lines.append(record.sequence)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _detect_sequence_length_outlier_rows(
    lengths_by_identifier: list[tuple[str, int]],
) -> list[SequenceLengthOutlier]:
    lengths = [length for _, length in lengths_by_identifier]
    if len(lengths) < 3:
        return []
    median_length = float(median(lengths))
    if median_length == 0.0:
        return []

    outliers: list[SequenceLengthOutlier] = []
    length_values = [float(length) for length in lengths]
    for identifier, raw_length in lengths_by_identifier:
        robust_z_score = _robust_z_score(raw_length, length_values)
        relative_deviation = abs(raw_length - median_length) / median_length
        if relative_deviation >= 0.2 or abs(robust_z_score or 0.0) >= 3.5:
            note = (
                "longer than baseline"
                if raw_length > median_length
                else "shorter than baseline"
            )
            outliers.append(
                SequenceLengthOutlier(
                    identifier=identifier,
                    raw_length=raw_length,
                    median_length=median_length,
                    robust_z_score=robust_z_score,
                    note=note,
                )
            )
    return sorted(outliers, key=lambda item: (item.raw_length, item.identifier))
