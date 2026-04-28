from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.core.alignment import AlignmentRecord, AlignmentSummary
from bijux_phylogenetics.errors import InvalidAlignmentError


def summarise_fasta(path: Path) -> AlignmentSummary:
    """Summarise a FASTA alignment without loading a heavy dependency."""
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
                    records.append(AlignmentRecord(identifier=current_identifier, sequence="".join(current_sequence)))
                current_identifier = line[1:].strip()
                current_sequence = []
                continue
            if current_identifier is None:
                raise InvalidAlignmentError(f"alignment sequence appears before any FASTA header in {path}")
            current_sequence.append(line)

    if current_identifier is not None:
        records.append(AlignmentRecord(identifier=current_identifier, sequence="".join(current_sequence)))

    if not records:
        raise InvalidAlignmentError(f"alignment contains no FASTA records: {path}")

    ids = [record.identifier for record in records]
    duplicate_ids = sorted(identifier for identifier in set(ids) if ids.count(identifier) > 1)
    if duplicate_ids:
        raise InvalidAlignmentError(f"alignment contains duplicate sequence ids: {', '.join(duplicate_ids)}")

    lengths = [len(record.sequence) for record in records]
    if min(lengths) != max(lengths):
        raise InvalidAlignmentError(
            f"alignment contains unequal sequence lengths: min={min(lengths)} max={max(lengths)}"
        )

    gap_characters = {"-"}
    missing_characters = {"?", "N", "n", "X", "x"}
    total_sites = len(records) * lengths[0]
    gap_count = sum(sum(1 for residue in record.sequence if residue in gap_characters) for record in records)
    missing_count = sum(sum(1 for residue in record.sequence if residue in missing_characters) for record in records)

    variable_site_count = 0
    parsimony_informative_site_count = 0
    for column in zip(*(record.sequence for record in records), strict=True):
        observed = [residue.upper() for residue in column if residue not in gap_characters and residue not in missing_characters]
        states = set(observed)
        if len(states) > 1:
            variable_site_count += 1
        if states and sum(observed.count(state) >= 2 for state in states) >= 2:
            parsimony_informative_site_count += 1

    return AlignmentSummary(
        path=path,
        sequence_count=len(records),
        alignment_length=lengths[0],
        min_sequence_length=min(lengths),
        max_sequence_length=max(lengths),
        ids=ids,
        missing_data_fraction=missing_count / total_sites,
        gap_fraction=gap_count / total_sites,
        variable_site_count=variable_site_count,
        parsimony_informative_site_count=parsimony_informative_site_count,
    )
