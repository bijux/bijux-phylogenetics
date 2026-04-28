from __future__ import annotations

from pathlib import Path
from statistics import median

from bijux_phylogenetics.core.alignment import (
    AlignmentAlphabet,
    AlignmentQualityReport,
    AlignmentLinkageReport,
    AlignmentRecord,
    AlignmentSummary,
    DuplicateSequenceGroup,
    InvalidAlignmentCharacter,
    NearDuplicateSequencePair,
    SequenceMissingness,
    SequenceCompositionOutlier,
    SequenceGCContent,
    SiteMissingness,
)
from bijux_phylogenetics.errors import AlignmentTaxonMismatchError, InvalidAlignmentError
from bijux_phylogenetics.io.trees import load_tree

_GAP_CHARACTERS = {"-"}
_MISSING_CHARACTERS = {"?", "N", "n", "X", "x"}
_DNA_CHARACTERS = set("ACGTNRYSWKMBDHVacgtnryswkmbdhv")
_RNA_CHARACTERS = set("ACGUNRYSWKMBDHVacgunryswkmbdhv")
_NUCLEOTIDE_GC_CHARACTERS = {"G", "C", "g", "c"}
_DNA_CANONICAL = ("A", "C", "G", "T")
_RNA_CANONICAL = ("A", "C", "G", "U")
_PROTEIN_CHARACTERS = set("ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz*")
_PROTEIN_CANONICAL = tuple("ACDEFGHIKLMNPQRSTVWY")


def _validate_fraction_threshold(threshold: float) -> None:
    if not 0.0 <= threshold <= 1.0:
        raise ValueError(f"threshold must be between 0 and 1 inclusive, got {threshold}")


def _observed_residues(records: list[AlignmentRecord]) -> list[str]:
    return [
        residue
        for record in records
        for residue in record.sequence
        if residue not in _GAP_CHARACTERS and residue not in _MISSING_CHARACTERS
    ]


def infer_alignment_alphabet(records: list[AlignmentRecord]) -> AlignmentAlphabet:
    """Infer whether an alignment is DNA, RNA, protein, or unknown."""
    observed = _observed_residues(records)
    if not observed:
        return "unknown"
    characters = set(observed)
    if characters <= _DNA_CHARACTERS and "U" not in {residue.upper() for residue in observed}:
        return "dna"
    if characters <= _RNA_CHARACTERS and "T" not in {residue.upper() for residue in observed}:
        return "rna"
    if characters <= _PROTEIN_CHARACTERS:
        return "protein"
    return "unknown"


def detect_invalid_alignment_characters(
    path: Path,
    *,
    alphabet: AlignmentAlphabet,
) -> list[InvalidAlignmentCharacter]:
    """List sequence characters invalid for the declared alphabet."""
    records = load_fasta_alignment(path)
    if alphabet == "dna":
        allowed = _DNA_CHARACTERS | _GAP_CHARACTERS | _MISSING_CHARACTERS
    elif alphabet == "rna":
        allowed = _RNA_CHARACTERS | _GAP_CHARACTERS | _MISSING_CHARACTERS
    elif alphabet == "protein":
        allowed = _PROTEIN_CHARACTERS | _GAP_CHARACTERS | _MISSING_CHARACTERS
    else:
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


def _normalized_frequency(values: list[str], alphabet: tuple[str, ...]) -> dict[str, float]:
    if not values:
        return {}
    total = len(values)
    return {
        character: round(sum(1 for value in values if value.upper() == character) / total, 15)
        for character in alphabet
        if any(value.upper() == character for value in values)
    }


def compute_nucleotide_composition(records: list[AlignmentRecord], *, alphabet: AlignmentAlphabet) -> dict[str, float]:
    """Compute canonical nucleotide composition for DNA or RNA alignments."""
    observed = _observed_residues(records)
    if alphabet == "dna":
        return _normalized_frequency(observed, _DNA_CANONICAL)
    if alphabet == "rna":
        return _normalized_frequency(observed, _RNA_CANONICAL)
    return {}


def compute_amino_acid_composition(records: list[AlignmentRecord], *, alphabet: AlignmentAlphabet) -> dict[str, float]:
    """Compute canonical amino-acid composition for protein alignments."""
    if alphabet != "protein":
        return {}
    return _normalized_frequency(_observed_residues(records), _PROTEIN_CANONICAL)


def _sequence_gc_fraction(sequence: str) -> float | None:
    comparable = [residue for residue in sequence if residue.upper() in {"A", "C", "G", "T", "U"}]
    if not comparable:
        return None
    return round(sum(1 for residue in comparable if residue in _NUCLEOTIDE_GC_CHARACTERS) / len(comparable), 15)


def compute_per_sequence_gc_content(records: list[AlignmentRecord], *, alphabet: AlignmentAlphabet) -> list[SequenceGCContent]:
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


def compute_whole_alignment_gc_content(records: list[AlignmentRecord], *, alphabet: AlignmentAlphabet) -> float | None:
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
    return round(sum(1 for residue in comparable if residue in _NUCLEOTIDE_GC_CHARACTERS) / len(comparable), 15)


def detect_composition_outlier_sequences(
    path: Path,
    *,
    deviation_threshold: float = 0.25,
) -> list[SequenceCompositionOutlier]:
    """Detect sequences with unusually deviant GC or amino-acid composition."""
    records = load_fasta_alignment(path)
    alphabet = infer_alignment_alphabet(records)
    if alphabet in {"dna", "rna"}:
        per_sequence_gc = [row for row in compute_per_sequence_gc_content(records, alphabet=alphabet) if row.gc_fraction is not None]
        if len(per_sequence_gc) < 2:
            return []
        baseline = median(row.gc_fraction for row in per_sequence_gc if row.gc_fraction is not None)
        return sorted(
            [
                SequenceCompositionOutlier(
                    identifier=row.identifier,
                    deviation=round(abs(float(row.gc_fraction) - baseline), 15),
                )
                for row in per_sequence_gc
                if row.gc_fraction is not None and abs(float(row.gc_fraction) - baseline) > deviation_threshold
            ],
            key=lambda item: (-item.deviation, item.identifier),
        )

    if alphabet == "protein":
        profile = compute_amino_acid_composition(records, alphabet=alphabet)
        outliers: list[SequenceCompositionOutlier] = []
        for record in records:
            sequence_profile = compute_amino_acid_composition([record], alphabet=alphabet)
            deviation = sum(abs(sequence_profile.get(key, 0.0) - profile.get(key, 0.0)) for key in set(sequence_profile) | set(profile))
            if deviation > deviation_threshold:
                outliers.append(
                    SequenceCompositionOutlier(identifier=record.identifier, deviation=round(deviation, 15))
                )
        return sorted(outliers, key=lambda item: (-item.deviation, item.identifier))
    return []


def detect_identical_duplicate_sequences(path: Path) -> list[DuplicateSequenceGroup]:
    """Group sequences that are exactly identical over the full aligned string."""
    records = load_fasta_alignment(path)
    grouped: dict[str, list[str]] = {}
    for record in records:
        grouped.setdefault(record.sequence, []).append(record.identifier)
    return [
        DuplicateSequenceGroup(identifiers=sorted(identifiers), sequence=sequence)
        for sequence, identifiers in sorted(grouped.items())
        if len(identifiers) > 1
    ]


def _pairwise_identity(left: str, right: str) -> tuple[float, int]:
    comparable_pairs = [
        (left_residue, right_residue)
        for left_residue, right_residue in zip(left, right, strict=True)
        if left_residue not in _GAP_CHARACTERS
        and right_residue not in _GAP_CHARACTERS
        and left_residue not in _MISSING_CHARACTERS
        and right_residue not in _MISSING_CHARACTERS
    ]
    if not comparable_pairs:
        return 0.0, 0
    matches = sum(1 for left_residue, right_residue in comparable_pairs if left_residue.upper() == right_residue.upper())
    comparable_sites = len(comparable_pairs)
    return round(matches / comparable_sites, 15), comparable_sites


def detect_near_duplicate_sequences(
    path: Path,
    *,
    identity_threshold: float,
) -> list[NearDuplicateSequencePair]:
    """Return sequence pairs above the given identity threshold."""
    _validate_fraction_threshold(identity_threshold)
    records = load_fasta_alignment(path)
    near_duplicates: list[NearDuplicateSequencePair] = []
    for index, left in enumerate(records):
        for right in records[index + 1 :]:
            identity, comparable_sites = _pairwise_identity(left.sequence, right.sequence)
            if comparable_sites > 0 and identity >= identity_threshold and left.sequence != right.sequence:
                near_duplicates.append(
                    NearDuplicateSequencePair(
                        left_identifier=left.identifier,
                        right_identifier=right.identifier,
                        identity=identity,
                        comparable_sites=comparable_sites,
                    )
                )
    return near_duplicates


def build_alignment_quality_report(path: Path) -> AlignmentQualityReport:
    """Generate a higher-level alignment quality report from composition and identity diagnostics."""
    summary = summarise_fasta(path)
    records = load_fasta_alignment(path)
    inferred_alphabet = summary.inferred_alphabet
    invalid_characters = (
        []
        if inferred_alphabet == "unknown"
        else detect_invalid_alignment_characters(path, alphabet=inferred_alphabet)
    )
    composition_outliers = detect_composition_outlier_sequences(path)
    duplicate_sequence_groups = detect_identical_duplicate_sequences(path)
    near_duplicate_pairs = detect_near_duplicate_sequences(path, identity_threshold=0.95)
    warnings: list[str] = []
    if invalid_characters:
        warnings.append("alignment contains characters invalid for the inferred alphabet")
    if composition_outliers:
        warnings.append("alignment contains composition outlier sequences")
    if duplicate_sequence_groups:
        warnings.append("alignment contains identical duplicate sequences")
    if near_duplicate_pairs:
        warnings.append("alignment contains near-duplicate sequences")
    return AlignmentQualityReport(
        path=path,
        sequence_count=summary.sequence_count,
        alignment_length=summary.alignment_length,
        missing_data_fraction=summary.missing_data_fraction,
        gap_fraction=summary.gap_fraction,
        variable_site_count=summary.variable_site_count,
        parsimony_informative_site_count=summary.parsimony_informative_site_count,
        inferred_alphabet=inferred_alphabet,
        invalid_characters=invalid_characters,
        composition_outliers=composition_outliers,
        duplicate_sequence_groups=duplicate_sequence_groups,
        near_duplicate_pairs=near_duplicate_pairs,
        warnings=warnings,
    )


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


def write_fasta_alignment(path: Path, records: list[AlignmentRecord]) -> Path:
    """Write FASTA records using a deterministic single-line sequence layout."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for record in records:
        lines.append(f">{record.identifier}")
        lines.append(record.sequence)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def summarise_fasta(path: Path) -> AlignmentSummary:
    """Summarise a FASTA alignment without loading a heavy dependency."""
    records = load_fasta_alignment(path)
    inferred_alphabet = infer_alignment_alphabet(records)
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

    invalid_characters = (
        []
        if inferred_alphabet == "unknown"
        else detect_invalid_alignment_characters(path, alphabet=inferred_alphabet)
    )
    nucleotide_composition = compute_nucleotide_composition(records, alphabet=inferred_alphabet)
    amino_acid_composition = compute_amino_acid_composition(records, alphabet=inferred_alphabet)
    per_sequence_gc_content = compute_per_sequence_gc_content(records, alphabet=inferred_alphabet)
    whole_alignment_gc_content = compute_whole_alignment_gc_content(records, alphabet=inferred_alphabet)
    composition_outliers = detect_composition_outlier_sequences(path)
    duplicate_sequence_groups = detect_identical_duplicate_sequences(path)
    near_duplicate_pairs = detect_near_duplicate_sequences(path, identity_threshold=0.95)

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
        inferred_alphabet=inferred_alphabet,
        invalid_characters=invalid_characters,
        nucleotide_composition=nucleotide_composition,
        amino_acid_composition=amino_acid_composition,
        per_sequence_gc_content=per_sequence_gc_content,
        whole_alignment_gc_content=whole_alignment_gc_content,
        composition_outliers=composition_outliers,
        duplicate_sequence_groups=duplicate_sequence_groups,
        near_duplicate_pairs=near_duplicate_pairs,
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


def detect_sequences_with_excessive_missing_data(
    path: Path,
    *,
    threshold: float,
) -> list[SequenceMissingness]:
    """Return sequences whose missing-data fraction exceeds the given threshold."""
    _validate_fraction_threshold(threshold)
    summary = summarise_fasta(path)
    return [row for row in summary.per_sequence_missingness if row.missing_fraction > threshold]


def detect_sites_with_excessive_missing_data(
    path: Path,
    *,
    threshold: float,
) -> list[SiteMissingness]:
    """Return alignment columns whose missing-data fraction exceeds the given threshold."""
    _validate_fraction_threshold(threshold)
    summary = summarise_fasta(path)
    return [row for row in summary.per_site_missingness if row.missing_fraction > threshold]
