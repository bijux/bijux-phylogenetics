from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    FastaDuplicateIdentifier,
    FastaEmptySequence,
    FastaIllegalCharacter,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from .character_policy import (
    _DNA_CHARACTERS,
    _EXPLICIT_MISSING_CHARACTERS,
    _GAP_CHARACTERS,
    _PROTEIN_CHARACTERS,
    _PROTEIN_EXCLUSIVE_CHARACTERS_UPPER,
    _RNA_CHARACTERS,
    _compatible_raw_sequence_types,
    _observed_raw_sequence_characters,
)


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
    if alphabet == "dna":
        allowed = _DNA_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
    elif alphabet == "rna":
        allowed = _RNA_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
    elif alphabet == "protein":
        allowed = _PROTEIN_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
    else:
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
