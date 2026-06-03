from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    InvalidAlignmentCharacter,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from .character_policy import _allowed_characters_for_sequence_type
from .raw_scan import _detect_duplicate_identifiers


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


def detect_invalid_alignment_characters(
    path: Path,
    *,
    alphabet: AlignmentAlphabet,
) -> list[InvalidAlignmentCharacter]:
    """List sequence characters invalid for the declared alphabet."""
    records = load_fasta_alignment(path)
    return _detect_invalid_alignment_characters_records(records, alphabet=alphabet)


def write_fasta_alignment(path: Path, records: list[AlignmentRecord]) -> Path:
    """Write FASTA records using a deterministic single-line sequence layout."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for record in records:
        lines.append(f">{record.identifier}")
        lines.append(record.sequence)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
