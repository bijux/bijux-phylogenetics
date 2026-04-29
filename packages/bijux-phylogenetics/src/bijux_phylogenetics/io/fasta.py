from __future__ import annotations

from pathlib import Path
from statistics import median

from Bio.Data import CodonTable

from bijux_phylogenetics.core.alignment import (
    AlignmentAlphabet,
    AlignmentMethodReadiness,
    AlignmentQualityReport,
    AlignmentReadinessReport,
    AlignmentLinkageReport,
    AlignmentRecord,
    AlignmentSequenceKindReport,
    AlignmentSummary,
    AlignmentSuspiciousRegion,
    AlignmentTrimReport,
    AlignmentWindowSummary,
    CodingAlignmentDiagnostics,
    DuplicateSequenceGroup,
    FrameshiftLikeSequence,
    InvalidAlignmentCharacter,
    NearDuplicateSequencePair,
    PairwiseSequenceIdentity,
    PartialCodonSequence,
    RemovedAlignmentSequence,
    SequenceMissingness,
    SequenceLengthOutlier,
    SequenceIdentityMatrix,
    SequenceCompositionOutlier,
    SequenceGCContent,
    SequenceUncertaintyProfile,
    SiteMissingness,
    SiteUncertaintyProfile,
    StopCodonObservation,
    TranslationReport,
    TrimmedAlignmentColumn,
)
from bijux_phylogenetics.errors import AlignmentTaxonMismatchError, InvalidAlignmentError
from bijux_phylogenetics.io.trees import load_tree

_GAP_CHARACTERS = {"-"}
_EXPLICIT_MISSING_CHARACTERS = {"?"}
_DNA_CHARACTERS = set("ACGTNRYSWKMBDHVacgtnryswkmbdhv")
_RNA_CHARACTERS = set("ACGUNRYSWKMBDHVacgunryswkmbdhv")
_NUCLEOTIDE_GC_CHARACTERS = {"G", "C", "g", "c"}
_DNA_CANONICAL = ("A", "C", "G", "T")
_RNA_CANONICAL = ("A", "C", "G", "U")
_PROTEIN_CHARACTERS = set("ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz*")
_PROTEIN_CANONICAL = tuple("ACDEFGHIKLMNPQRSTVWY")
_DNA_BASES = {"A", "C", "G", "T"}
_DNA_AMBIGUITY_UPPER = {"N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}
_RNA_AMBIGUITY_UPPER = {"N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}
_PROTEIN_AMBIGUITY_UPPER = {"B", "J", "X", "Z"}
_STANDARD_DNA_TABLE = CodonTable.unambiguous_dna_by_name["Standard"]
_STANDARD_FORWARD_TABLE = _STANDARD_DNA_TABLE.forward_table
_STANDARD_STOP_CODONS = set(_STANDARD_DNA_TABLE.stop_codons)


def _validate_fraction_threshold(threshold: float) -> None:
    if not 0.0 <= threshold <= 1.0:
        raise ValueError(f"threshold must be between 0 and 1 inclusive, got {threshold}")


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
    return _is_explicit_missing(residue) or _is_ambiguity_character(residue, alphabet=alphabet)


def _observed_residues(records: list[AlignmentRecord]) -> list[str]:
    return [
        residue
        for record in records
        for residue in record.sequence
        if residue not in _GAP_CHARACTERS and not _is_explicit_missing(residue)
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
        allowed = _DNA_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
    elif alphabet == "rna":
        allowed = _RNA_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
    elif alphabet == "protein":
        allowed = _PROTEIN_CHARACTERS | _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS
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
        gc_values = [float(row.gc_fraction) for row in per_sequence_gc if row.gc_fraction is not None]
        baseline = median(gc_values)
        return sorted(
            [
                SequenceCompositionOutlier(
                    identifier=row.identifier,
                    deviation=round(abs(float(row.gc_fraction) - baseline), 15),
                    robust_z_score=None if row.gc_fraction is None else _robust_z_score(float(row.gc_fraction), gc_values),
                )
                for row in per_sequence_gc
                if row.gc_fraction is not None
                and (
                    abs(float(row.gc_fraction) - baseline) > deviation_threshold
                    or abs(_robust_z_score(float(row.gc_fraction), gc_values) or 0.0) >= 3.5
                )
            ],
            key=lambda item: (-item.deviation, item.identifier),
        )

    if alphabet == "protein":
        profile = compute_amino_acid_composition(records, alphabet=alphabet)
        deviations_by_identifier: dict[str, float] = {}
        for record in records:
            sequence_profile = compute_amino_acid_composition([record], alphabet=alphabet)
            deviation = sum(abs(sequence_profile.get(key, 0.0) - profile.get(key, 0.0)) for key in set(sequence_profile) | set(profile))
            deviations_by_identifier[record.identifier] = round(deviation, 15)
        deviation_values = list(deviations_by_identifier.values())
        outliers: list[SequenceCompositionOutlier] = []
        for record in records:
            deviation = deviations_by_identifier[record.identifier]
            robust_z_score = _robust_z_score(deviation, deviation_values)
            if deviation > deviation_threshold or abs(robust_z_score or 0.0) >= 3.5:
                outliers.append(
                    SequenceCompositionOutlier(
                        identifier=record.identifier,
                        deviation=deviation,
                        robust_z_score=robust_z_score,
                    )
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
        and not _is_explicit_missing(left_residue)
        and not _is_explicit_missing(right_residue)
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


def _trim_columns(
    records: list[AlignmentRecord],
    *,
    keep_positions: list[int],
) -> list[AlignmentRecord]:
    return [
        AlignmentRecord(
            identifier=record.identifier,
            sequence="".join(record.sequence[index] for index in keep_positions),
        )
        for record in records
    ]


def remove_all_gap_columns(path: Path) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Remove columns composed entirely of gap characters."""
    records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    removed_positions = set(summary.all_gap_columns)
    keep_positions = [index for index in range(summary.alignment_length) if (index + 1) not in removed_positions]
    trimmed_records = _trim_columns(records, keep_positions=keep_positions)
    return trimmed_records, AlignmentTrimReport(
        path=path,
        original_sequence_count=summary.sequence_count,
        trimmed_sequence_count=len(trimmed_records),
        original_alignment_length=summary.alignment_length,
        trimmed_alignment_length=len(keep_positions),
        removed_columns=[
            TrimmedAlignmentColumn(position=position, reason="all-gap")
            for position in summary.all_gap_columns
        ],
        removed_sequences=[],
    )


def remove_all_missing_columns(path: Path) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Remove columns composed entirely of missing-data symbols."""
    records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    removed_positions = set(summary.all_missing_columns)
    keep_positions = [index for index in range(summary.alignment_length) if (index + 1) not in removed_positions]
    trimmed_records = _trim_columns(records, keep_positions=keep_positions)
    return trimmed_records, AlignmentTrimReport(
        path=path,
        original_sequence_count=summary.sequence_count,
        trimmed_sequence_count=len(trimmed_records),
        original_alignment_length=summary.alignment_length,
        trimmed_alignment_length=len(keep_positions),
        removed_columns=[
            TrimmedAlignmentColumn(position=position, reason="all-missing")
            for position in summary.all_missing_columns
        ],
        removed_sequences=[],
    )


def trim_columns_above_missingness_threshold(
    path: Path,
    *,
    threshold: float,
) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Remove columns whose missing-data fraction exceeds the given threshold."""
    _validate_fraction_threshold(threshold)
    records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    removed_positions = {
        row.position
        for row in detect_sites_with_excessive_missing_data(path, threshold=threshold)
    }
    keep_positions = [index for index in range(summary.alignment_length) if (index + 1) not in removed_positions]
    trimmed_records = _trim_columns(records, keep_positions=keep_positions)
    return trimmed_records, AlignmentTrimReport(
        path=path,
        original_sequence_count=summary.sequence_count,
        trimmed_sequence_count=len(trimmed_records),
        original_alignment_length=summary.alignment_length,
        trimmed_alignment_length=len(keep_positions),
        removed_columns=[
            TrimmedAlignmentColumn(position=position, reason="missingness-threshold")
            for position in sorted(removed_positions)
        ],
        removed_sequences=[],
    )


def remove_sequences_above_missingness_threshold(
    path: Path,
    *,
    threshold: float,
) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Remove sequences whose missing-data fraction exceeds the given threshold."""
    _validate_fraction_threshold(threshold)
    records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    excessive = {
        row.identifier: row.missing_fraction
        for row in detect_sequences_with_excessive_missing_data(path, threshold=threshold)
    }
    trimmed_records = [record for record in records if record.identifier not in excessive]
    removed_sequences = [
        RemovedAlignmentSequence(
            identifier=record.identifier,
            missing_fraction=excessive[record.identifier],
            reason="missingness-threshold",
        )
        for record in records
        if record.identifier in excessive
    ]
    return trimmed_records, AlignmentTrimReport(
        path=path,
        original_sequence_count=summary.sequence_count,
        trimmed_sequence_count=len(trimmed_records),
        original_alignment_length=summary.alignment_length,
        trimmed_alignment_length=summary.alignment_length,
        removed_columns=[],
        removed_sequences=removed_sequences,
    )


def trim_alignment(
    path: Path,
    *,
    remove_all_gap_sites: bool = True,
    remove_all_missing_sites: bool = True,
    site_missingness_threshold: float | None = None,
    sequence_missingness_threshold: float | None = None,
) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Trim an alignment with explicit deterministic transform reporting."""
    if site_missingness_threshold is not None:
        _validate_fraction_threshold(site_missingness_threshold)
    if sequence_missingness_threshold is not None:
        _validate_fraction_threshold(sequence_missingness_threshold)

    original_records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    records = list(original_records)
    removed_columns: list[TrimmedAlignmentColumn] = []
    removed_sequences: list[RemovedAlignmentSequence] = []

    if sequence_missingness_threshold is not None:
        excessive = {
            row.identifier: row.missing_fraction
            for row in detect_sequences_with_excessive_missing_data(
                path,
                threshold=sequence_missingness_threshold,
            )
        }
        removed_sequences.extend(
            RemovedAlignmentSequence(
                identifier=record.identifier,
                missing_fraction=excessive[record.identifier],
                reason="missingness-threshold",
            )
            for record in records
            if record.identifier in excessive
        )
        records = [record for record in records if record.identifier not in excessive]

    removed_positions: set[int] = set()
    if remove_all_gap_sites:
        removed_positions.update(summary.all_gap_columns)
        removed_columns.extend(
            TrimmedAlignmentColumn(position=position, reason="all-gap")
            for position in summary.all_gap_columns
        )
    if remove_all_missing_sites:
        removed_positions.update(
            position for position in summary.all_missing_columns if position not in removed_positions
        )
        removed_columns.extend(
            TrimmedAlignmentColumn(position=position, reason="all-missing")
            for position in summary.all_missing_columns
            if position not in summary.all_gap_columns or not remove_all_gap_sites
        )
    if site_missingness_threshold is not None:
        excessive_sites = [
            row.position
            for row in detect_sites_with_excessive_missing_data(
                path,
                threshold=site_missingness_threshold,
            )
            if row.position not in removed_positions
        ]
        removed_positions.update(excessive_sites)
        removed_columns.extend(
            TrimmedAlignmentColumn(position=position, reason="missingness-threshold")
            for position in excessive_sites
        )

    keep_positions = [index for index in range(summary.alignment_length) if (index + 1) not in removed_positions]
    trimmed_records = _trim_columns(records, keep_positions=keep_positions)
    removed_columns.sort(key=lambda item: item.position)
    removed_sequences.sort(key=lambda item: item.identifier)
    return trimmed_records, AlignmentTrimReport(
        path=path,
        original_sequence_count=summary.sequence_count,
        trimmed_sequence_count=len(trimmed_records),
        original_alignment_length=summary.alignment_length,
        trimmed_alignment_length=len(keep_positions),
        removed_columns=removed_columns,
        removed_sequences=removed_sequences,
    )


def compute_pairwise_sequence_identity_matrix(path: Path) -> SequenceIdentityMatrix:
    """Compute a deterministic pairwise sequence identity matrix."""
    records = load_fasta_alignment(path)
    pairs: list[PairwiseSequenceIdentity] = []
    for left_index, left in enumerate(records):
        for right_index, right in enumerate(records):
            if right_index < left_index:
                continue
            if left_index == right_index:
                pairs.append(
                    PairwiseSequenceIdentity(
                        left_identifier=left.identifier,
                        right_identifier=right.identifier,
                        identity=1.0,
                        comparable_sites=len(
                            [
                                residue
                                for residue in left.sequence
                                if residue not in _GAP_CHARACTERS and not _is_explicit_missing(residue)
                            ]
                        ),
                    )
                )
                continue
            identity, comparable_sites = _pairwise_identity(left.sequence, right.sequence)
            pairs.append(
                PairwiseSequenceIdentity(
                    left_identifier=left.identifier,
                    right_identifier=right.identifier,
                    identity=identity if comparable_sites > 0 else None,
                    comparable_sites=comparable_sites,
                )
            )
    return SequenceIdentityMatrix(
        path=path,
        identifiers=[record.identifier for record in records],
        pairs=pairs,
    )


def write_sequence_identity_matrix(path: Path, report: SequenceIdentityMatrix) -> Path:
    """Write a pairwise sequence identity matrix as a deterministic TSV."""
    rows = {(pair.left_identifier, pair.right_identifier): pair for pair in report.pairs}
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["left_identifier\tright_identifier\tidentity\tcomparable_sites"]
    for left in report.identifiers:
        for right in report.identifiers:
            pair = rows.get((left, right)) or rows.get((right, left))
            if pair is None:
                continue
            identity = "" if pair.identity is None else format(pair.identity, ".15g")
            lines.append(
                f"{left}\t{right}\t{identity}\t{pair.comparable_sites}"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _coding_residues(sequence: str) -> str:
    return "".join(
        residue.upper().replace("U", "T")
        for residue in sequence
        if residue not in _GAP_CHARACTERS and not _is_explicit_missing(residue)
    )


def detect_frameshift_like_sequences(path: Path) -> list[FrameshiftLikeSequence]:
    """Detect coding sequences whose comparable length is not divisible by three."""
    records = load_fasta_alignment(path)
    flagged: list[FrameshiftLikeSequence] = []
    for record in records:
        comparable_length = len(_coding_residues(record.sequence))
        remainder = comparable_length % 3
        if remainder != 0:
            flagged.append(
                FrameshiftLikeSequence(
                    identifier=record.identifier,
                    comparable_length=comparable_length,
                    remainder=remainder,
                )
            )
    return flagged


def _iter_codon_windows(sequence: str) -> list[tuple[int, str]]:
    return [
        (position, sequence[position - 1 : position + 2])
        for position in range(1, len(sequence) + 1, 3)
        if len(sequence[position - 1 : position + 2]) == 3
    ]


def detect_stop_codons(path: Path) -> list[StopCodonObservation]:
    """Detect stop codons in a coding alignment under the standard genetic code."""
    records = load_fasta_alignment(path)
    stop_codons: list[StopCodonObservation] = []
    for record in records:
        coding_sequence = _coding_residues(record.sequence)
        codons = _iter_codon_windows(coding_sequence)
        for codon_index, (start, codon) in enumerate(codons, start=1):
            if set(codon) <= _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS:
                continue
            if any(base not in _DNA_BASES for base in codon):
                continue
            if codon in _STANDARD_STOP_CODONS:
                stop_codons.append(
                    StopCodonObservation(
                        identifier=record.identifier,
                        codon_index=codon_index,
                        nucleotide_start=start,
                        codon=codon,
                        terminal=codon_index == len(codons),
                    )
                )
    return stop_codons


def inspect_coding_alignment(path: Path) -> CodingAlignmentDiagnostics:
    """Inspect one nucleotide alignment as a coding sequence dataset."""
    summary = summarise_fasta(path)
    if summary.inferred_alphabet not in {"dna", "rna"}:
        raise InvalidAlignmentError(
            f"coding diagnostics require a nucleotide alignment, got alphabet '{summary.inferred_alphabet}'"
        )
    partial_codon_sequences = [
        PartialCodonSequence(
            identifier=row.identifier,
            comparable_length=row.comparable_length,
            trailing_bases=row.remainder,
        )
        for row in detect_frameshift_like_sequences(path)
    ]
    return CodingAlignmentDiagnostics(
        path=path,
        sequence_count=summary.sequence_count,
        alignment_length=summary.alignment_length,
        alignment_length_multiple_of_three=summary.alignment_length % 3 == 0,
        frameshift_like_sequences=detect_frameshift_like_sequences(path),
        partial_codon_sequences=partial_codon_sequences,
        stop_codons=detect_stop_codons(path),
    )


def _translate_codon(codon: str) -> str:
    normalized = codon.upper().replace("U", "T")
    if set(normalized) <= _GAP_CHARACTERS:
        return "-"
    if any(base in _GAP_CHARACTERS or _is_explicit_missing(base) for base in normalized):
        return "X"
    if any(base not in _DNA_BASES for base in normalized):
        return "X"
    if normalized in _STANDARD_STOP_CODONS:
        return "*"
    return _STANDARD_FORWARD_TABLE.get(normalized, "X")


def translate_coding_alignment(path: Path) -> tuple[list[AlignmentRecord], TranslationReport]:
    """Translate an aligned nucleotide coding sequence dataset to amino acids."""
    summary = summarise_fasta(path)
    if summary.inferred_alphabet not in {"dna", "rna"}:
        raise InvalidAlignmentError(
            f"coding translation requires a nucleotide alignment, got alphabet '{summary.inferred_alphabet}'"
        )
    if summary.alignment_length % 3 != 0:
        raise InvalidAlignmentError(
            f"alignment length must be divisible by 3 for coding translation, got {summary.alignment_length}"
        )

    records = load_fasta_alignment(path)
    translated_records = [
        AlignmentRecord(
            identifier=record.identifier,
            sequence="".join(_translate_codon(codon) for _, codon in _iter_codon_windows(record.sequence)),
        )
        for record in records
    ]
    diagnostics = inspect_coding_alignment(path)
    return translated_records, TranslationReport(
        source_path=path,
        translated_sequence_count=len(translated_records),
        source_alignment_length=summary.alignment_length,
        translated_alignment_length=summary.alignment_length // 3,
        stop_codon_count=len(diagnostics.stop_codons),
        frameshift_like_sequence_count=len(diagnostics.frameshift_like_sequences),
    )


def summarize_alignment_windows(
    path: Path,
    *,
    window_size: int = 30,
    step_size: int = 10,
) -> list[AlignmentWindowSummary]:
    """Summarize an aligned FASTA file in sliding windows."""
    if window_size <= 0:
        raise ValueError(f"window_size must be positive, got {window_size}")
    if step_size <= 0:
        raise ValueError(f"step_size must be positive, got {step_size}")

    summary = summarise_fasta(path)
    records = load_fasta_alignment(path)
    windows: list[AlignmentWindowSummary] = []
    for start_index in range(0, summary.alignment_length, step_size):
        end_index = min(start_index + window_size, summary.alignment_length)
        if end_index <= start_index:
            continue
        # Flatten into per-position columns for metrics.
        columns = list(zip(*(record.sequence[start_index:end_index] for record in records), strict=True))
        if not columns:
            continue

        gap_count = 0
        missing_count = 0
        ambiguity_count = 0
        variable_sites = 0
        disagreement_fractions: list[float] = []
        comparable_sites = 0
        total_residues = len(records) * len(columns)
        alphabet = summary.inferred_alphabet
        for column in columns:
            gap_count += sum(1 for residue in column if residue in _GAP_CHARACTERS)
            missing_count += sum(1 for residue in column if _is_explicit_missing(residue))
            ambiguity_count += sum(
                1 for residue in column if _is_ambiguity_character(residue, alphabet=alphabet)
            )
            comparable = [
                residue.upper()
                for residue in column
                if residue not in _GAP_CHARACTERS
                and not _is_explicit_missing(residue)
                and not _is_ambiguity_character(residue, alphabet=alphabet)
            ]
            if comparable:
                comparable_sites += len(comparable)
                state_counts = {state: comparable.count(state) for state in set(comparable)}
                if len(state_counts) > 1:
                    variable_sites += 1
                disagreement_fractions.append(1.0 - (max(state_counts.values()) / len(comparable)))
            else:
                disagreement_fractions.append(0.0)

        windows.append(
            AlignmentWindowSummary(
                start=start_index + 1,
                end=end_index,
                site_count=len(columns),
                gap_fraction=round(gap_count / total_residues, 15),
                missing_fraction=round(missing_count / total_residues, 15),
                ambiguity_fraction=round(ambiguity_count / total_residues, 15),
                variable_fraction=round(variable_sites / len(columns), 15),
                disagreement_fraction=round(sum(disagreement_fractions) / len(disagreement_fractions), 15),
                comparable_fraction=round(comparable_sites / total_residues, 15),
            )
        )
        if end_index == summary.alignment_length:
            break
    return windows


def detect_over_aligned_regions(
    path: Path,
    *,
    window_size: int = 30,
    step_size: int = 10,
) -> list[AlignmentSuspiciousRegion]:
    """Flag suspicious windows that look excessively gap-heavy or over-regularized."""
    regions: list[AlignmentSuspiciousRegion] = []
    for window in summarize_alignment_windows(path, window_size=window_size, step_size=step_size):
        uncertainty_fraction = window.gap_fraction + window.missing_fraction + window.ambiguity_fraction
        if uncertainty_fraction >= 0.4 and window.variable_fraction <= 0.2:
            regions.append(
                AlignmentSuspiciousRegion(
                    start=window.start,
                    end=window.end,
                    kind="over_aligned",
                    score=round(uncertainty_fraction - window.variable_fraction, 15),
                    note="gap- or ambiguity-heavy window with little residual variability; review for aggressive or artifactual alignment",
                )
            )
    return regions


def detect_under_aligned_regions(
    path: Path,
    *,
    window_size: int = 30,
    step_size: int = 10,
) -> list[AlignmentSuspiciousRegion]:
    """Flag suspicious windows with strong local mismatch or gap disorder."""
    regions: list[AlignmentSuspiciousRegion] = []
    for window in summarize_alignment_windows(path, window_size=window_size, step_size=step_size):
        if window.variable_fraction >= 0.7 and window.disagreement_fraction >= 0.35:
            regions.append(
                AlignmentSuspiciousRegion(
                    start=window.start,
                    end=window.end,
                    kind="under_aligned",
                    score=round(window.variable_fraction + window.disagreement_fraction, 15),
                    note="high local mismatch and disagreement suggest the region may require realignment or masking",
                )
            )
    return regions


def summarize_alignment_readiness(path: Path) -> AlignmentReadinessReport:
    """Classify whether an input alignment is ready for key downstream analysis families."""
    sequence_kind = classify_alignment_sequences(path)
    records = load_fasta_records(path)
    inferred_alphabet = infer_alignment_alphabet(records)
    sequence_count = len(records)
    alignment_length = sequence_kind.max_sequence_length if sequence_kind.state != "raw_sequence_fasta" else None
    methods: list[AlignmentMethodReadiness] = []
    warnings: list[str] = []
    length_outliers = detect_sequence_length_outliers(path)

    if sequence_kind.state == "raw_sequence_fasta":
        common_blocker = ["input sequences are not yet aligned"]
        methods.extend(
            [
                AlignmentMethodReadiness(analysis="distance", ready=False, blockers=common_blocker, warnings=[]),
                AlignmentMethodReadiness(analysis="maximum_likelihood", ready=False, blockers=common_blocker, warnings=[]),
                AlignmentMethodReadiness(analysis="bayesian", ready=False, blockers=common_blocker, warnings=[]),
                AlignmentMethodReadiness(analysis="coding", ready=False, blockers=common_blocker, warnings=[]),
                AlignmentMethodReadiness(analysis="protein", ready=inferred_alphabet == "protein", blockers=[] if inferred_alphabet == "protein" else common_blocker, warnings=[]),
            ]
        )
        if length_outliers:
            warnings.append("raw sequences include substantial length outliers before alignment")
        return AlignmentReadinessReport(
            path=path,
            sequence_kind=sequence_kind,
            inferred_alphabet=inferred_alphabet,
            sequence_count=sequence_count,
            alignment_length=alignment_length,
            methods=methods,
            warnings=warnings,
        )

    summary = summarise_fasta(path)
    quality = build_alignment_quality_report(path)
    over_aligned = detect_over_aligned_regions(path)
    under_aligned = detect_under_aligned_regions(path)
    coding = inspect_coding_alignment(path) if inferred_alphabet in {"dna", "rna"} else None

    def _method(analysis: str, ready: bool, blockers: list[str], extra_warnings: list[str]) -> AlignmentMethodReadiness:
        return AlignmentMethodReadiness(analysis=analysis, ready=ready, blockers=blockers, warnings=extra_warnings)

    generic_warnings: list[str] = []
    if sequence_kind.state == "ambiguous_equal_length_fasta":
        generic_warnings.append("equal-length ungapped FASTA may be aligned, but prior alignment cannot be proved from sequence shape alone")
    if over_aligned:
        generic_warnings.append("one or more windows look suspiciously over-aligned")
    if under_aligned:
        generic_warnings.append("one or more windows look suspiciously under-aligned")
    if quality.composition_outliers:
        generic_warnings.append("composition outliers may bias downstream inference")
    if quality.sequence_length_outliers:
        generic_warnings.append("sequence length outliers suggest problematic trimming, concatenation, or mixed loci")
    warnings.extend(generic_warnings)

    base_alignment_blockers = ["alignment contains invalid characters for the inferred alphabet"] if quality.invalid_characters else []
    distance_blockers = list(base_alignment_blockers)
    if summary.variable_site_count == 0:
        distance_blockers.append("alignment has no variable sites")
    methods.append(
        _method(
            "distance",
            ready=not distance_blockers,
            blockers=distance_blockers,
            extra_warnings=generic_warnings,
        )
    )

    model_blockers = list(base_alignment_blockers)
    if summary.variable_site_count == 0:
        model_blockers.append("alignment has no variable sites")
    methods.append(
        _method(
            "maximum_likelihood",
            ready=not model_blockers,
            blockers=model_blockers,
            extra_warnings=generic_warnings,
        )
    )
    methods.append(
        _method(
            "bayesian",
            ready=not model_blockers,
            blockers=model_blockers,
            extra_warnings=generic_warnings,
        )
    )

    coding_blockers: list[str] = []
    coding_warnings = list(generic_warnings)
    if inferred_alphabet not in {"dna", "rna"}:
        coding_blockers.append("coding analysis requires a nucleotide alignment")
    elif coding is not None:
        if not coding.alignment_length_multiple_of_three:
            coding_blockers.append("alignment length is not divisible by three")
        if coding.frameshift_like_sequences:
            coding_blockers.append("one or more sequences contain partial codons after gaps and missing data are removed")
        if any(not stop.terminal for stop in coding.stop_codons):
            coding_blockers.append("one or more sequences contain premature stop codons")
        if any(stop.terminal for stop in coding.stop_codons):
            coding_warnings.append("terminal stop codons were detected and should be verified against coding conventions")
    methods.append(_method("coding", ready=not coding_blockers, blockers=coding_blockers, extra_warnings=coding_warnings))

    protein_blockers = [] if inferred_alphabet == "protein" else ["protein analysis requires an amino-acid alignment"]
    methods.append(_method("protein", ready=not protein_blockers, blockers=protein_blockers, extra_warnings=generic_warnings))

    return AlignmentReadinessReport(
        path=path,
        sequence_kind=sequence_kind,
        inferred_alphabet=inferred_alphabet,
        sequence_count=summary.sequence_count,
        alignment_length=summary.alignment_length,
        methods=methods,
        warnings=warnings,
    )


def build_alignment_quality_report(path: Path) -> AlignmentQualityReport:
    """Generate a higher-level alignment quality report from composition and identity diagnostics."""
    summary = summarise_fasta(path)
    inferred_alphabet = summary.inferred_alphabet
    invalid_characters = (
        []
        if inferred_alphabet == "unknown"
        else detect_invalid_alignment_characters(path, alphabet=inferred_alphabet)
    )
    composition_outliers = detect_composition_outlier_sequences(path)
    sequence_length_outliers = detect_sequence_length_outliers(path)
    duplicate_sequence_groups = detect_identical_duplicate_sequences(path)
    near_duplicate_pairs = detect_near_duplicate_sequences(path, identity_threshold=0.95)
    warnings: list[str] = []
    if invalid_characters:
        warnings.append("alignment contains characters invalid for the inferred alphabet")
    if composition_outliers:
        warnings.append("alignment contains composition outlier sequences")
    if sequence_length_outliers:
        warnings.append("alignment contains raw-sequence length outliers")
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
        ambiguity_fraction=summary.ambiguity_fraction,
        variable_site_count=summary.variable_site_count,
        parsimony_informative_site_count=summary.parsimony_informative_site_count,
        inferred_alphabet=inferred_alphabet,
        invalid_characters=invalid_characters,
        composition_outliers=composition_outliers,
        sequence_length_outliers=sequence_length_outliers,
        duplicate_sequence_groups=duplicate_sequence_groups,
        near_duplicate_pairs=near_duplicate_pairs,
        warnings=warnings,
    )


def load_fasta_records(path: Path) -> list[AlignmentRecord]:
    """Load FASTA records without imposing equal-length alignment requirements."""
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

    return records


def classify_alignment_sequences(path: Path) -> AlignmentSequenceKindReport:
    """Classify whether a FASTA input already behaves like an aligned dataset."""
    records = load_fasta_records(path)
    lengths = [len(record.sequence) for record in records]
    min_length = min(lengths)
    max_length = max(lengths)
    has_gaps = any(any(residue in _GAP_CHARACTERS for residue in record.sequence) for record in records)

    if min_length != max_length:
        state = "raw_sequence_fasta"
        note = "sequence lengths differ, so the FASTA cannot yet be treated as one aligned matrix"
    elif has_gaps:
        state = "aligned"
        note = "equal sequence lengths with gap characters are consistent with an existing alignment"
    else:
        state = "ambiguous_equal_length_fasta"
        note = "sequence lengths are equal but no gap characters were observed, so prior alignment cannot be proved"

    return AlignmentSequenceKindReport(
        path=path,
        sequence_count=len(records),
        min_sequence_length=min_length,
        max_sequence_length=max_length,
        has_gap_characters=has_gaps,
        state=state,
        note=note,
    )


def detect_sequence_length_outliers(path: Path) -> list[SequenceLengthOutlier]:
    """Detect unusually short or long raw sequences before alignment assumptions are imposed."""
    records = load_fasta_records(path)
    lengths = [len(record.sequence) for record in records]
    if len(lengths) < 3:
        return []
    median_length = float(median(lengths))
    if median_length == 0.0:
        return []

    outliers: list[SequenceLengthOutlier] = []
    for record in records:
        robust_z_score = _robust_z_score(len(record.sequence), [float(length) for length in lengths])
        relative_deviation = abs(len(record.sequence) - median_length) / median_length
        if relative_deviation >= 0.2 or abs(robust_z_score or 0.0) >= 3.5:
            note = "longer than baseline" if len(record.sequence) > median_length else "shorter than baseline"
            outliers.append(
                SequenceLengthOutlier(
                    identifier=record.identifier,
                    raw_length=len(record.sequence),
                    median_length=median_length,
                    robust_z_score=robust_z_score,
                    note=note,
                )
            )
    return sorted(outliers, key=lambda item: (item.raw_length, item.identifier))


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


def summarise_fasta(path: Path) -> AlignmentSummary:
    """Summarise a FASTA alignment without loading a heavy dependency."""
    records = load_fasta_alignment(path)
    inferred_alphabet = infer_alignment_alphabet(records)
    ids = [record.identifier for record in records]
    lengths = [len(record.sequence) for record in records]
    total_sites = len(records) * lengths[0]
    gap_count = sum(sum(1 for residue in record.sequence if residue in _GAP_CHARACTERS) for record in records)
    missing_count = sum(
        sum(1 for residue in record.sequence if _is_missing_like(residue, alphabet=inferred_alphabet))
        for record in records
    )
    ambiguity_count = sum(
        sum(1 for residue in record.sequence if _is_ambiguity_character(residue, alphabet=inferred_alphabet))
        for record in records
    )
    per_sequence_missingness = [
        SequenceMissingness(
            identifier=record.identifier,
            missing_fraction=sum(
                1 for residue in record.sequence if _is_missing_like(residue, alphabet=inferred_alphabet)
            )
            / lengths[0],
        )
        for record in records
    ]
    per_sequence_uncertainty = [
        SequenceUncertaintyProfile(
            identifier=record.identifier,
            gap_fraction=sum(1 for residue in record.sequence if residue in _GAP_CHARACTERS) / lengths[0],
            missing_fraction=sum(
                1 for residue in record.sequence if _is_missing_like(residue, alphabet=inferred_alphabet)
            )
            / lengths[0],
            ambiguity_fraction=sum(
                1 for residue in record.sequence if _is_ambiguity_character(residue, alphabet=inferred_alphabet)
            )
            / lengths[0],
        )
        for record in records
    ]
    per_site_missingness: list[SiteMissingness] = []
    per_site_uncertainty: list[SiteUncertaintyProfile] = []
    all_gap_columns: list[int] = []
    all_missing_columns: list[int] = []

    constant_site_count = 0
    variable_site_count = 0
    parsimony_informative_site_count = 0
    for position, column in enumerate(zip(*(record.sequence for record in records), strict=True), start=1):
        missing_fraction = sum(
            1 for residue in column if _is_missing_like(residue, alphabet=inferred_alphabet)
        ) / len(records)
        per_site_missingness.append(SiteMissingness(position=position, missing_fraction=missing_fraction))
        per_site_uncertainty.append(
            SiteUncertaintyProfile(
                position=position,
                gap_fraction=sum(1 for residue in column if residue in _GAP_CHARACTERS) / len(records),
                missing_fraction=missing_fraction,
                ambiguity_fraction=sum(
                    1 for residue in column if _is_ambiguity_character(residue, alphabet=inferred_alphabet)
                )
                / len(records),
            )
        )
        if all(residue in _GAP_CHARACTERS for residue in column):
            all_gap_columns.append(position)
        if all(_is_missing_like(residue, alphabet=inferred_alphabet) for residue in column):
            all_missing_columns.append(position)
        observed = [
            residue.upper()
            for residue in column
            if residue not in _GAP_CHARACTERS
            and not _is_explicit_missing(residue)
            and not _is_ambiguity_character(residue, alphabet=inferred_alphabet)
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
        ambiguity_fraction=ambiguity_count / total_sites,
        per_sequence_missingness=per_sequence_missingness,
        per_sequence_uncertainty=per_sequence_uncertainty,
        per_site_missingness=per_site_missingness,
        per_site_uncertainty=per_site_uncertainty,
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
