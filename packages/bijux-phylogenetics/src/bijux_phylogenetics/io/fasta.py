from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.core.alignment import (
    AlignmentLinkageReport,
    AlignmentRecord,
    AlignmentSummary,
    SequenceMissingness,
    SiteMissingness,
)
from bijux_phylogenetics.errors import AlignmentTaxonMismatchError, InvalidAlignmentError
from bijux_phylogenetics.io.trees import load_tree

_GAP_CHARACTERS = {"-"}
_MISSING_CHARACTERS = {"?", "N", "n", "X", "x"}


def load_fasta_alignment(path: Path) -> list[AlignmentRecord]:
    """Load FASTA records using the repository's deterministic alignment contract."""
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

    return records


def summarise_fasta(path: Path) -> AlignmentSummary:
    """Summarise a FASTA alignment without loading a heavy dependency."""
    records = load_fasta_alignment(path)
    ids = [record.identifier for record in records]
    lengths = [len(record.sequence) for record in records]
    total_sites = len(records) * lengths[0]
    gap_count = sum(sum(1 for residue in record.sequence if residue in _GAP_CHARACTERS) for record in records)
    missing_count = sum(sum(1 for residue in record.sequence if residue in _MISSING_CHARACTERS) for record in records)
    per_sequence_missingness = [
        SequenceMissingness(
            identifier=record.identifier,
            missing_fraction=sum(1 for residue in record.sequence if residue in _MISSING_CHARACTERS) / lengths[0],
        )
        for record in records
    ]
    per_site_missingness: list[SiteMissingness] = []
    all_gap_columns: list[int] = []
    all_missing_columns: list[int] = []

    constant_site_count = 0
    variable_site_count = 0
    parsimony_informative_site_count = 0
    for position, column in enumerate(zip(*(record.sequence for record in records), strict=True), start=1):
        missing_fraction = sum(1 for residue in column if residue in _MISSING_CHARACTERS) / len(records)
        per_site_missingness.append(SiteMissingness(position=position, missing_fraction=missing_fraction))
        if all(residue in _GAP_CHARACTERS for residue in column):
            all_gap_columns.append(position)
        if all(residue in _MISSING_CHARACTERS for residue in column):
            all_missing_columns.append(position)
        observed = [
            residue.upper()
            for residue in column
            if residue not in _GAP_CHARACTERS and residue not in _MISSING_CHARACTERS
        ]
        states = set(observed)
        if len(states) == 1 and observed:
            constant_site_count += 1
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
        per_sequence_missingness=per_sequence_missingness,
        per_site_missingness=per_site_missingness,
        all_gap_columns=all_gap_columns,
        all_missing_columns=all_missing_columns,
        constant_site_count=constant_site_count,
        variable_site_count=variable_site_count,
        parsimony_informative_site_count=parsimony_informative_site_count,
    )


def link_alignment_to_tree(
    tree_path: Path,
    alignment_path: Path,
    *,
    strict: bool = False,
) -> AlignmentLinkageReport:
    """Report how alignment sequence identifiers join against a tree."""
    tree = load_tree(tree_path)
    alignment = summarise_fasta(alignment_path)
    tree_taxa = set(tree.tip_names)
    alignment_ids = set(alignment.ids)
    missing_from_alignment = sorted(tree_taxa - alignment_ids)
    extra_alignment_ids = sorted(alignment_ids - tree_taxa)

    if strict and (missing_from_alignment or extra_alignment_ids):
        raise AlignmentTaxonMismatchError(
            "alignment linkage mismatch: "
            f"{len(missing_from_alignment)} tree taxa missing from alignment and "
            f"{len(extra_alignment_ids)} alignment ids absent from tree"
        )

    usable_taxa = sorted(tree_taxa & alignment_ids)
    return AlignmentLinkageReport(
        tree_path=tree_path,
        alignment_path=alignment_path,
        tree_taxa=len(tree_taxa),
        alignment_ids=len(alignment_ids),
        linked_taxa=len(usable_taxa),
        usable_taxa=usable_taxa,
        missing_from_alignment=missing_from_alignment,
        extra_alignment_ids=extra_alignment_ids,
    )
