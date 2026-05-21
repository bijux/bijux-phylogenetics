from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    InvalidAlignmentCharacter,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from .character_policy import (
    _allowed_characters_for_sequence_type as _allowed_characters_for_sequence_type,
    _compatible_raw_sequence_types as _compatible_raw_sequence_types,
    _DNA_CANONICAL as _DNA_CANONICAL,
    _DNA_CHARACTERS as _DNA_CHARACTERS,
    _DNA_CHARACTERS_UPPER as _DNA_CHARACTERS_UPPER,
    _EXPLICIT_MISSING_CHARACTERS as _EXPLICIT_MISSING_CHARACTERS,
    _GAP_CHARACTERS as _GAP_CHARACTERS,
    _is_ambiguity_character as _is_ambiguity_character,
    _is_explicit_missing as _is_explicit_missing,
    _is_missing_like as _is_missing_like,
    _NUCLEOTIDE_GC_CHARACTERS as _NUCLEOTIDE_GC_CHARACTERS,
    _observed_raw_sequence_characters as _observed_raw_sequence_characters,
    _observed_residues as _observed_residues,
    _ordered_sequence_types as _ordered_sequence_types,
    _PROTEIN_CANONICAL as _PROTEIN_CANONICAL,
    _PROTEIN_CHARACTERS as _PROTEIN_CHARACTERS,
    _PROTEIN_CHARACTERS_UPPER as _PROTEIN_CHARACTERS_UPPER,
    _PROTEIN_EXCLUSIVE_CHARACTERS_UPPER as _PROTEIN_EXCLUSIVE_CHARACTERS_UPPER,
    _RNA_CANONICAL as _RNA_CANONICAL,
    _RNA_CHARACTERS as _RNA_CHARACTERS,
    _RNA_CHARACTERS_UPPER as _RNA_CHARACTERS_UPPER,
    _SUPPORTED_SEQUENCE_CHARACTERS_UPPER as _SUPPORTED_SEQUENCE_CHARACTERS_UPPER,
    _validate_fraction_threshold as _validate_fraction_threshold,
)
from .raw_scan import (
    _detect_duplicate_identifiers as _detect_duplicate_identifiers,
    _detect_empty_sequences as _detect_empty_sequences,
    _detect_illegal_sequence_characters as _detect_illegal_sequence_characters,
    _RawFastaScan as _RawFastaScan,
    _scan_raw_fasta as _scan_raw_fasta,
)
from .input_validation import (
    _build_fasta_input_summary_from_scan as _build_fasta_input_summary_from_scan,
    _build_fasta_sequence_type_report as _build_fasta_sequence_type_report,
    _detect_sequence_length_outlier_rows as _detect_sequence_length_outlier_rows,
    detect_fasta_sequence_type as detect_fasta_sequence_type,
    _median_absolute_deviation as _median_absolute_deviation,
    _robust_z_score as _robust_z_score,
)
from .composition import (
    compute_amino_acid_composition as compute_amino_acid_composition,
    compute_nucleotide_composition as compute_nucleotide_composition,
    compute_per_sequence_gc_content as compute_per_sequence_gc_content,
    compute_whole_alignment_gc_content as compute_whole_alignment_gc_content,
    infer_alignment_alphabet as infer_alignment_alphabet,
    _normalized_frequency as _normalized_frequency,
    _sequence_gc_fraction as _sequence_gc_fraction,
)


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
