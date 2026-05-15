from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from statistics import median

from Bio.Data import CodonTable

from bijux_phylogenetics.core.alignment import (
    AlignmentBaseFrequencyReport,
    AlignmentAlphabet,
    AlignmentAmbiguousColumnReport,
    AlignmentCleaningReport,
    AlignmentComparisonReport,
    AlignmentCompositionShift,
    AlignmentFilterProfile,
    AlignmentForensicReport,
    AlignmentGroupRetention,
    AlignmentLinkageReport,
    AlignmentLowInformationReport,
    AlignmentMethodReadiness,
    AlignmentMissingDataConcentration,
    AlignmentQualityReport,
    AlignmentReadinessReport,
    AlignmentRecord,
    AlignmentSequenceKindReport,
    AlignmentSignalWarning,
    AlignmentSummary,
    AlignmentSuspiciousRegion,
    AlignmentTrimReport,
    AlignmentWindowSummary,
    AmbiguousAlignmentColumn,
    CodingAlignmentDiagnostics,
    CodingSequenceExclusion,
    CodingSequencePreparationReport,
    DuplicateSequenceGroup,
    DuplicateSequencePolicyAction,
    DuplicateSequencePolicyReport,
    FastaDuplicateIdentifier,
    FastaEmptySequence,
    FastaIdentifierRepair,
    FastaIllegalCharacter,
    FastaInputSummary,
    FastaInputValidationReport,
    FastaRemovedRecord,
    FastaRepairReport,
    FastaSequenceTypeReport,
    FrameshiftLikeSequence,
    InvalidAlignmentCharacter,
    InvalidCodonObservation,
    NearDuplicateSequencePair,
    NucleotideStateFrequencyRow,
    AlignmentSegregatingSiteReport,
    PairwiseSequenceIdentity,
    PartialCodonSequence,
    RemovedAlignmentSequence,
    SegregatingSiteRow,
    SequenceCodingBehavior,
    SequenceCompositionOutlier,
    SequenceGCContent,
    SequenceIdentityMatrix,
    SequenceLengthOutlier,
    SequenceMissingness,
    SequenceQualityRankingReport,
    SequenceQualityRankingRow,
    SequenceUncertaintyProfile,
    SiteMissingness,
    SiteUncertaintyProfile,
    StopCodonObservation,
    TranslationCodonObservation,
    TranslationReport,
    TranslationSequenceExclusion,
    TrimmedAlignmentColumn,
)
from bijux_phylogenetics.core.metadata import load_taxon_table
from bijux_phylogenetics.errors import (
    AlignmentTaxonMismatchError,
    InvalidAlignmentError,
)
from bijux_phylogenetics.io.trees import load_tree

_GAP_CHARACTERS = {"-"}
_EXPLICIT_MISSING_CHARACTERS = {"?"}
_DNA_CHARACTERS = set("ACGTNRYSWKMBDHVacgtnryswkmbdhv")
_RNA_CHARACTERS = set("ACGUNRYSWKMBDHVacgunryswkmbdhv")
_NUCLEOTIDE_GC_CHARACTERS = {"G", "C", "g", "c"}
_DNA_CANONICAL = ("A", "C", "G", "T")
_RNA_CANONICAL = ("A", "C", "G", "U")
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
_PROTEIN_CANONICAL = tuple("ACDEFGHIKLMNPQRSTVWY")
_DNA_BASES = {"A", "C", "G", "T"}
_DNA_AMBIGUITY_UPPER = {"N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}
_RNA_AMBIGUITY_UPPER = {"N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}
_PROTEIN_AMBIGUITY_UPPER = {"B", "J", "X", "Z"}
_APE_DNA_STATE_ORDER = (
    "a",
    "c",
    "g",
    "t",
    "r",
    "m",
    "w",
    "s",
    "k",
    "y",
    "v",
    "h",
    "d",
    "b",
    "n",
    "-",
    "?",
)
_APE_SEGREGATING_STATE_SETS: dict[str, frozenset[str]] = {
    "A": frozenset({"A"}),
    "C": frozenset({"C"}),
    "G": frozenset({"G"}),
    "T": frozenset({"T"}),
    "M": frozenset({"A", "C"}),
    "R": frozenset({"A", "G"}),
    "W": frozenset({"A", "T"}),
    "S": frozenset({"C", "G"}),
    "K": frozenset({"G", "T"}),
    "Y": frozenset({"C", "T"}),
    "V": frozenset({"A", "C", "G"}),
    "H": frozenset({"A", "C", "T"}),
    "D": frozenset({"A", "G", "T"}),
    "B": frozenset({"C", "G", "T"}),
    "N": frozenset({"A", "C", "G", "T"}),
}
_DEFAULT_GENETIC_CODE_ID = 1
_LOW_INFORMATION_SITE_THRESHOLD = 1
_LOW_INFORMATION_FRACTION_THRESHOLD = 0.0
_MAX_DEFAULT_NEAR_DUPLICATE_SEQUENCE_COUNT = 256
_ALIGNMENT_FILTER_PROFILES: dict[str, AlignmentFilterProfile] = {
    "conservative": AlignmentFilterProfile(
        name="conservative",
        remove_all_gap_sites=True,
        remove_all_missing_sites=True,
        site_missingness_threshold=0.8,
        sequence_missingness_threshold=0.6,
        preserve_codon_structure=False,
        note="retain most data while removing only the most uncertainty-heavy sites and sequences",
    ),
    "moderate": AlignmentFilterProfile(
        name="moderate",
        remove_all_gap_sites=True,
        remove_all_missing_sites=True,
        site_missingness_threshold=0.5,
        sequence_missingness_threshold=0.4,
        preserve_codon_structure=False,
        note="balanced trimming for routine phylogenetic alignment cleaning",
    ),
    "aggressive": AlignmentFilterProfile(
        name="aggressive",
        remove_all_gap_sites=True,
        remove_all_missing_sites=True,
        site_missingness_threshold=0.3,
        sequence_missingness_threshold=0.25,
        preserve_codon_structure=False,
        note="remove strongly uncertainty-heavy sites and weak sequences even at the cost of data loss",
    ),
    "coding-safe": AlignmentFilterProfile(
        name="coding-safe",
        remove_all_gap_sites=True,
        remove_all_missing_sites=True,
        site_missingness_threshold=0.5,
        sequence_missingness_threshold=0.35,
        preserve_codon_structure=True,
        note="trim codon-aligned nucleotide data while preserving codon phase by removing full codons together",
    ),
    "phylogenomics-scale": AlignmentFilterProfile(
        name="phylogenomics-scale",
        remove_all_gap_sites=True,
        remove_all_missing_sites=True,
        site_missingness_threshold=0.7,
        sequence_missingness_threshold=0.7,
        preserve_codon_structure=False,
        note="favor taxon retention for large concatenated matrices while removing empty or highly degraded regions",
    ),
}


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


def _resolve_genetic_code_table(
    genetic_code: int | str | None,
) -> tuple[int, str, dict[str, str], set[str]]:
    if genetic_code is None:
        code_id = _DEFAULT_GENETIC_CODE_ID
        table = CodonTable.unambiguous_dna_by_id[code_id]
        return code_id, table.names[0], table.forward_table, set(table.stop_codons)
    if isinstance(genetic_code, int):
        code_id = genetic_code
        try:
            table = CodonTable.unambiguous_dna_by_id[code_id]
        except KeyError as error:
            raise InvalidAlignmentError(
                f"unsupported genetic code id for coding analysis: {code_id}"
            ) from error
        return code_id, table.names[0], table.forward_table, set(table.stop_codons)
    text = genetic_code.strip()
    if not text:
        raise InvalidAlignmentError("genetic code name must not be empty")
    if text.isdigit():
        return _resolve_genetic_code_table(int(text))
    lowered = text.lower()
    for table in CodonTable.unambiguous_dna_by_id.values():
        if lowered in {name.lower() for name in table.names}:
            return (
                int(table.id),
                table.names[0],
                table.forward_table,
                set(table.stop_codons),
            )
    raise InvalidAlignmentError(
        f"unsupported genetic code for coding analysis: {genetic_code}"
    )


def _validate_fraction_threshold(threshold: float) -> None:
    if not 0.0 <= threshold <= 1.0:
        raise ValueError(
            f"threshold must be between 0 and 1 inclusive, got {threshold}"
        )


def list_alignment_filter_profiles() -> list[AlignmentFilterProfile]:
    """Return the supported named alignment-cleaning profiles."""
    return list(_ALIGNMENT_FILTER_PROFILES.values())


def get_alignment_filter_profile(name: str) -> AlignmentFilterProfile:
    """Resolve one named alignment-cleaning profile."""
    try:
        return _ALIGNMENT_FILTER_PROFILES[name]
    except KeyError as error:
        available = ", ".join(sorted(_ALIGNMENT_FILTER_PROFILES))
        raise ValueError(
            f"unknown alignment filter profile '{name}', expected one of: {available}"
        ) from error


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
        residue, alphabet=alphabet
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


def detect_invalid_alignment_characters(
    path: Path,
    *,
    alphabet: AlignmentAlphabet,
) -> list[InvalidAlignmentCharacter]:
    """List sequence characters invalid for the declared alphabet."""
    records = load_fasta_alignment(path)
    return _detect_invalid_alignment_characters_records(records, alphabet=alphabet)


def _detect_invalid_alignment_characters_records(
    records: list[AlignmentRecord],
    *,
    alphabet: AlignmentAlphabet,
) -> list[InvalidAlignmentCharacter]:
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


def _normalized_frequency(
    values: list[str], alphabet: tuple[str, ...]
) -> dict[str, float]:
    if not values:
        return {}
    total = len(values)
    return {
        character: round(
            sum(1 for value in values if value.upper() == character) / total, 15
        )
        for character in alphabet
        if any(value.upper() == character for value in values)
    }


def compute_nucleotide_composition(
    records: list[AlignmentRecord], *, alphabet: AlignmentAlphabet
) -> dict[str, float]:
    """Compute canonical nucleotide composition for DNA or RNA alignments."""
    observed = _observed_residues(records)
    if alphabet == "dna":
        return _normalized_frequency(observed, _DNA_CANONICAL)
    if alphabet == "rna":
        return _normalized_frequency(observed, _RNA_CANONICAL)
    return {}


def compute_amino_acid_composition(
    records: list[AlignmentRecord], *, alphabet: AlignmentAlphabet
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
    records: list[AlignmentRecord], *, alphabet: AlignmentAlphabet
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
    records: list[AlignmentRecord], *, alphabet: AlignmentAlphabet
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


def _normalize_ape_nucleotide_state(residue: str) -> str | None:
    normalized = residue.lower().replace("u", "t")
    if normalized in _APE_DNA_STATE_ORDER:
        return normalized
    return None


def _ape_nucleotide_state_counts(sequence: str) -> dict[str, int]:
    counts = {state: 0 for state in _APE_DNA_STATE_ORDER}
    for residue in sequence:
        normalized = _normalize_ape_nucleotide_state(residue)
        if normalized is None:
            continue
        counts[normalized] += 1
    return counts


def compute_alignment_base_frequency_report(path: Path) -> AlignmentBaseFrequencyReport:
    """Compute ape-style nucleotide state frequencies for one DNA or RNA alignment."""
    records = load_fasta_alignment(path)
    alphabet = infer_alignment_alphabet(records)
    if alphabet not in {"dna", "rna"}:
        if any(
            _normalize_ape_nucleotide_state(residue) is None
            for record in records
            for residue in record.sequence
        ):
            raise InvalidAlignmentError(
                "ape-style nucleotide base frequencies require a dna or rna alignment"
            )
        alphabet = "dna"

    alignment_counts = {state: 0 for state in _APE_DNA_STATE_ORDER}
    per_sequence_rows: list[NucleotideStateFrequencyRow] = []
    for record in records:
        sequence_counts = _ape_nucleotide_state_counts(record.sequence)
        sequence_total = sum(sequence_counts.values())
        for state in _APE_DNA_STATE_ORDER:
            alignment_counts[state] += sequence_counts[state]
            per_sequence_rows.append(
                NucleotideStateFrequencyRow(
                    scope="sequence",
                    identifier=record.identifier,
                    state=state,
                    count=sequence_counts[state],
                    frequency=(
                        0.0
                        if sequence_total == 0
                        else round(sequence_counts[state] / sequence_total, 15)
                    ),
                )
            )

    alignment_total = sum(alignment_counts.values())
    alignment_rows = [
        NucleotideStateFrequencyRow(
            scope="alignment",
            identifier=None,
            state=state,
            count=alignment_counts[state],
            frequency=(
                0.0
                if alignment_total == 0
                else round(alignment_counts[state] / alignment_total, 15)
            ),
        )
        for state in _APE_DNA_STATE_ORDER
    ]
    warnings: list[str] = []
    canonical_total = sum(
        alignment_counts[state] for state in _APE_DNA_STATE_ORDER[:4]
    )
    if canonical_total == 0:
        warnings.append(
            "alignment contains no canonical A/C/G/T residues, so ape-style base frequencies reflect only ambiguity, gap, and missing states"
        )
    return AlignmentBaseFrequencyReport(
        path=path,
        inferred_alphabet=alphabet,
        sequence_count=len(records),
        alignment_length=0 if not records else len(records[0].sequence),
        ambiguity_policy="count ambiguity codes as literal states",
        gap_policy="count gap characters as literal states",
        missing_data_policy="count explicit missing characters as literal states",
        state_order=list(_APE_DNA_STATE_ORDER),
        alignment_rows=alignment_rows,
        per_sequence_rows=per_sequence_rows,
        composition_outliers=_detect_composition_outlier_sequences_records(records),
        warnings=warnings,
    )


def write_alignment_base_frequency_table(
    path: Path,
    report: AlignmentBaseFrequencyReport,
) -> Path:
    """Write ape-style alignment and per-sequence nucleotide state frequencies as TSV."""
    rows = report.alignment_rows + report.per_sequence_rows
    lines = ["scope\tidentifier\tstate\tcount\tfrequency"]
    lines.extend(
        "\t".join(
            [
                row.scope,
                "" if row.identifier is None else row.identifier,
                row.state,
                str(row.count),
                format(row.frequency, ".15g"),
            ]
        )
        for row in rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _normalize_ape_segregating_state(residue: str) -> str | None:
    normalized = residue.upper().replace("U", "T")
    if normalized in _APE_SEGREGATING_STATE_SETS or normalized in {"-", "?"}:
        return normalized
    return None


def _leading_trailing_gaps_to_n(sequence: str) -> str:
    characters = list(sequence)
    left = 0
    while left < len(characters) and characters[left] == "-":
        characters[left] = "N"
        left += 1
    right = len(characters) - 1
    while right >= 0 and characters[right] == "-":
        characters[right] = "N"
        right -= 1
    return "".join(characters)


def _ape_segregating_states_different(left: str, right: str) -> bool:
    if left == "?" or right == "?":
        return False
    left_states = _APE_SEGREGATING_STATE_SETS.get(left)
    right_states = _APE_SEGREGATING_STATE_SETS.get(right)
    if left_states is None or right_states is None:
        return False
    return left_states.isdisjoint(right_states)


def _is_ape_known_base(state: str) -> bool:
    return state in {"A", "C", "G", "T"}


def _is_ape_gap_state(state: str) -> bool:
    return state == "-"


def _is_ape_missing_state(state: str) -> bool:
    return state == "?"


def _is_ape_segregating_column(column: list[str]) -> bool:
    if len(column) <= 1:
        return False
    index = 0
    end = len(column) - 1
    base = column[index]

    while not _is_ape_known_base(base):
        index += 1
        if index > end:
            return False
        current = column[index]
        if base != current:
            if not _is_ape_missing_state(base) and not _is_ape_missing_state(current):
                if not _is_ape_gap_state(base):
                    if _is_ape_gap_state(current):
                        return True
                    if _ape_segregating_states_different(current, base):
                        return True
                else:
                    return True
            base = current

    index += 1
    while index <= end:
        current = column[index]
        if current != base:
            if _is_ape_gap_state(current):
                return True
            if _ape_segregating_states_different(current, base):
                return True
        index += 1
    return False


def compute_alignment_segregating_site_report(
    path: Path,
) -> AlignmentSegregatingSiteReport:
    """Compute ape-style segregating sites for one DNA or RNA alignment."""
    records = load_fasta_alignment(path)
    alphabet = infer_alignment_alphabet(records)
    if alphabet not in {"dna", "rna"}:
        if any(
            _normalize_ape_segregating_state(residue) is None
            for record in records
            for residue in record.sequence
        ):
            raise InvalidAlignmentError(
                "ape-style segregating-site detection requires a dna or rna alignment"
            )
        alphabet = "dna"

    original_sequences = [
        "".join(
            _normalize_ape_segregating_state(residue) or ""
            for residue in record.sequence
        )
        for record in records
    ]
    effective_sequences = [
        _leading_trailing_gaps_to_n(sequence) for sequence in original_sequences
    ]

    segregating_site_positions: list[int] = []
    rows: list[SegregatingSiteRow] = []
    for position, (original_column, effective_column) in enumerate(
        zip(
            zip(*original_sequences, strict=True),
            zip(*effective_sequences, strict=True),
            strict=True,
        ),
        start=1,
    ):
        effective_states = list(effective_column)
        if not _is_ape_segregating_column(effective_states):
            continue
        segregating_site_positions.append(position)
        rows.append(
            SegregatingSiteRow(
                position=position,
                original_states="|".join(original_column),
                effective_states="|".join(effective_states),
                known_state_count=sum(
                    1 for state in effective_states if _is_ape_known_base(state)
                ),
                ambiguity_state_count=sum(
                    1
                    for state in effective_states
                    if state in _APE_SEGREGATING_STATE_SETS
                    and not _is_ape_known_base(state)
                ),
                gap_count=sum(1 for state in effective_states if _is_ape_gap_state(state)),
                missing_count=sum(
                    1 for state in effective_states if _is_ape_missing_state(state)
                ),
            )
        )

    warnings: list[str] = []
    canonical_total = sum(
        1
        for sequence in effective_sequences
        for state in sequence
        if _is_ape_known_base(state)
    )
    if canonical_total == 0:
        warnings.append(
            "alignment contains no canonical A/C/G/T residues, so ape-style segregating-site detection can only reflect ambiguity, gap, and missing states"
        )

    return AlignmentSegregatingSiteReport(
        path=path,
        inferred_alphabet=alphabet,
        sequence_count=len(records),
        alignment_length=0 if not records else len(records[0].sequence),
        ambiguity_policy="ambiguity states segregate only when they are surely incompatible with another observed state",
        gap_policy="internal gap characters can create segregating sites against known or incompatible ambiguous states",
        missing_data_policy="explicit missing characters do not create segregating sites",
        trailing_gap_policy="leading and trailing gap runs are normalized to N before ape-style segregating-site detection",
        segregating_site_positions=segregating_site_positions,
        rows=rows,
        warnings=warnings,
    )


def write_alignment_segregating_site_table(
    path: Path,
    report: AlignmentSegregatingSiteReport,
) -> Path:
    """Write one reviewer-facing ape-style segregating-site ledger as TSV."""
    lines = [
        "position\toriginal_states\teffective_states\tknown_state_count\tambiguity_state_count\tgap_count\tmissing_count"
    ]
    lines.extend(
        "\t".join(
            [
                str(row.position),
                row.original_states,
                row.effective_states,
                str(row.known_state_count),
                str(row.ambiguity_state_count),
                str(row.gap_count),
                str(row.missing_count),
            ]
        )
        for row in report.rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def detect_composition_outlier_sequences(
    path: Path,
    *,
    deviation_threshold: float = 0.25,
) -> list[SequenceCompositionOutlier]:
    """Detect sequences with unusually deviant GC or amino-acid composition."""
    records = load_fasta_alignment(path)
    return _detect_composition_outlier_sequences_records(
        records,
        deviation_threshold=deviation_threshold,
    )


def _detect_composition_outlier_sequences_records(
    records: list[AlignmentRecord],
    *,
    deviation_threshold: float = 0.25,
) -> list[SequenceCompositionOutlier]:
    alphabet = infer_alignment_alphabet(records)
    if alphabet in {"dna", "rna"}:
        per_sequence_gc = [
            row
            for row in compute_per_sequence_gc_content(records, alphabet=alphabet)
            if row.gc_fraction is not None
        ]
        if len(per_sequence_gc) < 2:
            return []
        gc_values = [
            float(row.gc_fraction)
            for row in per_sequence_gc
            if row.gc_fraction is not None
        ]
        baseline = median(gc_values)
        return sorted(
            [
                SequenceCompositionOutlier(
                    identifier=row.identifier,
                    deviation=round(abs(float(row.gc_fraction) - baseline), 15),
                    robust_z_score=None
                    if row.gc_fraction is None
                    else _robust_z_score(float(row.gc_fraction), gc_values),
                )
                for row in per_sequence_gc
                if row.gc_fraction is not None
                and (
                    abs(float(row.gc_fraction) - baseline) > deviation_threshold
                    or abs(_robust_z_score(float(row.gc_fraction), gc_values) or 0.0)
                    >= 3.5
                )
            ],
            key=lambda item: (-item.deviation, item.identifier),
        )

    if alphabet == "protein":
        profile = compute_amino_acid_composition(records, alphabet=alphabet)
        deviations_by_identifier: dict[str, float] = {}
        for record in records:
            sequence_profile = compute_amino_acid_composition(
                [record], alphabet=alphabet
            )
            deviation = sum(
                abs(sequence_profile.get(key, 0.0) - profile.get(key, 0.0))
                for key in set(sequence_profile) | set(profile)
            )
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
    return _detect_identical_duplicate_sequences_records(records)


def _detect_identical_duplicate_sequences_records(
    records: list[AlignmentRecord],
) -> list[DuplicateSequenceGroup]:
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
    matches = sum(
        1
        for left_residue, right_residue in comparable_pairs
        if left_residue.upper() == right_residue.upper()
    )
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
    return _detect_near_duplicate_sequences_records(
        records,
        identity_threshold=identity_threshold,
    )


def _detect_near_duplicate_sequences_records(
    records: list[AlignmentRecord],
    *,
    identity_threshold: float,
) -> list[NearDuplicateSequencePair]:
    near_duplicates: list[NearDuplicateSequencePair] = []
    for index, left in enumerate(records):
        for right in records[index + 1 :]:
            identity, comparable_sites = _pairwise_identity(
                left.sequence, right.sequence
            )
            if (
                comparable_sites > 0
                and identity >= identity_threshold
                and left.sequence != right.sequence
            ):
                near_duplicates.append(
                    NearDuplicateSequencePair(
                        left_identifier=left.identifier,
                        right_identifier=right.identifier,
                        identity=identity,
                        comparable_sites=comparable_sites,
                    )
                )
    return near_duplicates


def _collect_near_duplicate_sequences_for_summary(
    records: list[AlignmentRecord],
    *,
    identity_threshold: float,
) -> tuple[list[NearDuplicateSequencePair], bool]:
    if len(records) > _MAX_DEFAULT_NEAR_DUPLICATE_SEQUENCE_COUNT:
        return [], False
    return (
        _detect_near_duplicate_sequences_records(
            records,
            identity_threshold=identity_threshold,
        ),
        True,
    )


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


def _expand_removed_positions_to_groups(
    removed_positions: set[int],
    *,
    alignment_length: int,
    group_size: int,
) -> set[int]:
    if group_size <= 1:
        return set(removed_positions)
    expanded: set[int] = set()
    for position in removed_positions:
        group_start = ((position - 1) // group_size) * group_size + 1
        for grouped_position in range(
            group_start, min(group_start + group_size, alignment_length + 1)
        ):
            expanded.add(grouped_position)
    return expanded


def remove_all_gap_columns(
    path: Path,
) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Remove columns composed entirely of gap characters."""
    records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    removed_positions = set(summary.all_gap_columns)
    keep_positions = [
        index
        for index in range(summary.alignment_length)
        if (index + 1) not in removed_positions
    ]
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


def remove_all_missing_columns(
    path: Path,
) -> tuple[list[AlignmentRecord], AlignmentTrimReport]:
    """Remove columns composed entirely of missing-data symbols."""
    records = load_fasta_alignment(path)
    summary = summarise_fasta(path)
    removed_positions = set(summary.all_missing_columns)
    keep_positions = [
        index
        for index in range(summary.alignment_length)
        if (index + 1) not in removed_positions
    ]
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
    keep_positions = [
        index
        for index in range(summary.alignment_length)
        if (index + 1) not in removed_positions
    ]
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
        for row in detect_sequences_with_excessive_missing_data(
            path, threshold=threshold
        )
    }
    trimmed_records = [
        record for record in records if record.identifier not in excessive
    ]
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
    preserve_codon_structure: bool = False,
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
            position
            for position in summary.all_missing_columns
            if position not in removed_positions
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

    if preserve_codon_structure and summary.inferred_alphabet in {"dna", "rna"}:
        expanded_positions = _expand_removed_positions_to_groups(
            removed_positions,
            alignment_length=summary.alignment_length,
            group_size=3,
        )
        codon_phase_positions = sorted(expanded_positions - removed_positions)
        removed_positions = expanded_positions
        removed_columns.extend(
            TrimmedAlignmentColumn(position=position, reason="codon-phase-preservation")
            for position in codon_phase_positions
        )

    keep_positions = [
        index
        for index in range(summary.alignment_length)
        if (index + 1) not in removed_positions
    ]
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


def _composition_for_comparison(summary: AlignmentSummary) -> dict[str, float]:
    if summary.nucleotide_composition:
        return summary.nucleotide_composition
    return summary.amino_acid_composition


def compare_alignment_versions(
    left_path: Path, right_path: Path
) -> AlignmentComparisonReport:
    """Compare two alignments over taxa, sites, uncertainty burden, and composition."""
    left = summarise_fasta(left_path)
    right = summarise_fasta(right_path)
    composition_shifts: list[AlignmentCompositionShift] = []
    left_composition = _composition_for_comparison(left)
    right_composition = _composition_for_comparison(right)
    for component in sorted(set(left_composition) | set(right_composition)):
        before = left_composition.get(component, 0.0)
        after = right_composition.get(component, 0.0)
        composition_shifts.append(
            AlignmentCompositionShift(
                component=component,
                before=before,
                after=after,
                delta=round(after - before, 15),
            )
        )

    warnings: list[str] = []
    if left.inferred_alphabet != right.inferred_alphabet:
        warnings.append("alignment versions infer different alphabets")
    if set(left.ids) != set(right.ids):
        warnings.append("alignment versions contain different taxon sets")

    return AlignmentComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=sorted(set(left.ids) & set(right.ids)),
        left_only_taxa=sorted(set(left.ids) - set(right.ids)),
        right_only_taxa=sorted(set(right.ids) - set(left.ids)),
        left_alignment_length=left.alignment_length,
        right_alignment_length=right.alignment_length,
        left_missing_data_fraction=left.missing_data_fraction,
        right_missing_data_fraction=right.missing_data_fraction,
        left_gap_fraction=left.gap_fraction,
        right_gap_fraction=right.gap_fraction,
        left_variable_site_count=left.variable_site_count,
        right_variable_site_count=right.variable_site_count,
        left_parsimony_informative_site_count=left.parsimony_informative_site_count,
        right_parsimony_informative_site_count=right.parsimony_informative_site_count,
        composition_shifts=composition_shifts,
        warnings=warnings,
    )


def compare_alignment_summaries(
    left_path: Path,
    left: AlignmentSummary,
    right: AlignmentSummary,
) -> AlignmentComparisonReport:
    """Compare two already-computed alignment summaries."""
    composition_shifts: list[AlignmentCompositionShift] = []
    left_composition = _composition_for_comparison(left)
    right_composition = _composition_for_comparison(right)
    for component in sorted(set(left_composition) | set(right_composition)):
        before = left_composition.get(component, 0.0)
        after = right_composition.get(component, 0.0)
        composition_shifts.append(
            AlignmentCompositionShift(
                component=component,
                before=before,
                after=after,
                delta=round(after - before, 15),
            )
        )
    warnings: list[str] = []
    if left.inferred_alphabet != right.inferred_alphabet:
        warnings.append("alignment versions infer different alphabets")
    if set(left.ids) != set(right.ids):
        warnings.append("alignment versions contain different taxon sets")
    return AlignmentComparisonReport(
        left_path=left_path,
        right_path=right.path,
        shared_taxa=sorted(set(left.ids) & set(right.ids)),
        left_only_taxa=sorted(set(left.ids) - set(right.ids)),
        right_only_taxa=sorted(set(right.ids) - set(left.ids)),
        left_alignment_length=left.alignment_length,
        right_alignment_length=right.alignment_length,
        left_missing_data_fraction=left.missing_data_fraction,
        right_missing_data_fraction=right.missing_data_fraction,
        left_gap_fraction=left.gap_fraction,
        right_gap_fraction=right.gap_fraction,
        left_variable_site_count=left.variable_site_count,
        right_variable_site_count=right.variable_site_count,
        left_parsimony_informative_site_count=left.parsimony_informative_site_count,
        right_parsimony_informative_site_count=right.parsimony_informative_site_count,
        composition_shifts=composition_shifts,
        warnings=warnings,
    )


def _signal_warnings_for_cleaning(
    original: AlignmentSummary,
    cleaned: AlignmentSummary,
) -> list[AlignmentSignalWarning]:
    warnings: list[AlignmentSignalWarning] = []
    if original.parsimony_informative_site_count > 0:
        retained_fraction = (
            cleaned.parsimony_informative_site_count
            / original.parsimony_informative_site_count
        )
        if retained_fraction < 0.5:
            warnings.append(
                AlignmentSignalWarning(
                    code="informative-site-loss",
                    message="cleaning retained less than half of the original parsimony-informative sites",
                )
            )
    if cleaned.sequence_count < max(2, original.sequence_count // 2):
        warnings.append(
            AlignmentSignalWarning(
                code="taxon-coverage-loss",
                message="cleaning removed more than half of the original taxa",
            )
        )
    if cleaned.variable_site_count == 0 and original.variable_site_count > 0:
        warnings.append(
            AlignmentSignalWarning(
                code="variable-site-collapse",
                message="cleaning removed all variable sites from the alignment",
            )
        )
    return warnings


def _group_retention_after_cleaning(
    original_taxa: list[str],
    retained_taxa: list[str],
    *,
    table_path: Path | None,
    group_columns: list[str] | None,
) -> list[AlignmentGroupRetention]:
    if table_path is None or not group_columns:
        return []
    table = load_taxon_table(table_path)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    original_taxa_set = set(original_taxa)
    retained_taxa_set = set(retained_taxa)
    reports: list[AlignmentGroupRetention] = []
    for column in group_columns:
        if column not in table.columns:
            raise ValueError(f"group column '{column}' is not present in {table_path}")
        values = sorted(
            {
                rows_by_taxon[taxon][column]
                for taxon in original_taxa
                if taxon in rows_by_taxon and rows_by_taxon[taxon][column]
            }
        )
        for value in values:
            original_count = sum(
                1
                for taxon in original_taxa_set
                if taxon in rows_by_taxon and rows_by_taxon[taxon][column] == value
            )
            retained_count = sum(
                1
                for taxon in retained_taxa_set
                if taxon in rows_by_taxon and rows_by_taxon[taxon][column] == value
            )
            removed_count = original_count - retained_count
            reports.append(
                AlignmentGroupRetention(
                    column=column,
                    value=value,
                    original_count=original_count,
                    retained_count=retained_count,
                    removed_count=removed_count,
                    removed_fraction=0.0
                    if original_count == 0
                    else round(removed_count / original_count, 15),
                )
            )
    return reports


def clean_alignment_with_profile(
    path: Path,
    *,
    profile_name: str,
    group_table_path: Path | None = None,
    group_columns: list[str] | None = None,
) -> tuple[list[AlignmentRecord], AlignmentCleaningReport]:
    """Clean an alignment using one named profile and compare the result against the original."""
    profile = get_alignment_filter_profile(profile_name)
    cleaned_records, trim_report = trim_alignment(
        path,
        remove_all_gap_sites=profile.remove_all_gap_sites,
        remove_all_missing_sites=profile.remove_all_missing_sites,
        site_missingness_threshold=profile.site_missingness_threshold,
        sequence_missingness_threshold=profile.sequence_missingness_threshold,
        preserve_codon_structure=profile.preserve_codon_structure,
    )

    original = summarise_fasta(path)
    cleaned_summary = summarise_records_as_alignment_summary(
        path=path, records=cleaned_records
    )
    comparison = compare_alignment_summaries(path, original, cleaned_summary)
    signal_warnings = _signal_warnings_for_cleaning(original, cleaned_summary)
    group_retention = _group_retention_after_cleaning(
        original.ids,
        cleaned_summary.ids,
        table_path=group_table_path,
        group_columns=group_columns,
    )
    warnings = [warning.message for warning in signal_warnings]
    if any(
        group.removed_fraction >= 0.5 and group.original_count >= 1
        for group in group_retention
    ):
        warnings.append(
            "cleaning removed most taxa from one or more metadata or trait groups"
        )
    return cleaned_records, AlignmentCleaningReport(
        profile=profile,
        trim=trim_report,
        comparison=comparison,
        signal_warnings=signal_warnings,
        group_retention=group_retention,
        warnings=warnings,
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
                                if residue not in _GAP_CHARACTERS
                                and not _is_explicit_missing(residue)
                            ]
                        ),
                    )
                )
                continue
            identity, comparable_sites = _pairwise_identity(
                left.sequence, right.sequence
            )
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
    rows = {
        (pair.left_identifier, pair.right_identifier): pair for pair in report.pairs
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["left_identifier\tright_identifier\tidentity\tcomparable_sites"]
    for left in report.identifiers:
        for right in report.identifiers:
            pair = rows.get((left, right)) or rows.get((right, left))
            if pair is None:
                continue
            identity = "" if pair.identity is None else format(pair.identity, ".15g")
            lines.append(f"{left}\t{right}\t{identity}\t{pair.comparable_sites}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _coding_residues(sequence: str) -> str:
    return "".join(
        residue.upper().replace("U", "T")
        for residue in sequence
        if residue not in _GAP_CHARACTERS and not _is_explicit_missing(residue)
    )


def _detect_frameshift_like_records(
    records: list[AlignmentRecord],
) -> list[FrameshiftLikeSequence]:
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


def detect_frameshift_like_sequences(path: Path) -> list[FrameshiftLikeSequence]:
    """Detect coding sequences whose comparable length is not divisible by three."""
    records = load_fasta_alignment(path)
    return _detect_frameshift_like_records(records)


def _iter_codon_windows(sequence: str) -> list[tuple[int, str]]:
    return [
        (position, sequence[position - 1 : position + 2])
        for position in range(1, len(sequence) + 1, 3)
        if len(sequence[position - 1 : position + 2]) == 3
    ]


def _invalid_codon_reason(codon: str) -> str | None:
    normalized = codon.upper().replace("U", "T")
    if set(normalized) <= _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS:
        return None
    if any(
        base in _GAP_CHARACTERS or _is_explicit_missing(base) for base in normalized
    ):
        return "partial-missing-codon"
    unsupported = sorted(
        {base for base in normalized if base not in _DNA_CHARACTERS_UPPER}
    )
    if unsupported:
        return f"unsupported-residue:{''.join(unsupported)}"
    if any(base not in _DNA_BASES for base in normalized):
        return "ambiguous-codon"
    return None


def _detect_invalid_codons_in_records(
    records: list[AlignmentRecord],
) -> list[InvalidCodonObservation]:
    invalid_codons: list[InvalidCodonObservation] = []
    for record in records:
        coding_sequence = _coding_residues(record.sequence)
        for codon_index, (start, codon) in enumerate(
            _iter_codon_windows(coding_sequence),
            start=1,
        ):
            reason = _invalid_codon_reason(codon)
            if reason is None:
                continue
            invalid_codons.append(
                InvalidCodonObservation(
                    identifier=record.identifier,
                    codon_index=codon_index,
                    nucleotide_start=start,
                    codon=codon.upper(),
                    reason=reason,
                )
            )
    return invalid_codons


def _detect_stop_codons_in_records(
    records: list[AlignmentRecord],
    *,
    genetic_code: int | str | None = None,
) -> list[StopCodonObservation]:
    _code_id, _code_name, _forward_table, stop_codon_set = _resolve_genetic_code_table(
        genetic_code
    )
    stop_observations: list[StopCodonObservation] = []
    for record in records:
        coding_sequence = _coding_residues(record.sequence)
        codons = _iter_codon_windows(coding_sequence)
        for codon_index, (start, codon) in enumerate(codons, start=1):
            normalized = codon.upper().replace("U", "T")
            if _invalid_codon_reason(normalized) is not None:
                continue
            if normalized in stop_codon_set:
                stop_observations.append(
                    StopCodonObservation(
                        identifier=record.identifier,
                        codon_index=codon_index,
                        nucleotide_start=start,
                        codon=normalized,
                        terminal=codon_index == len(codons),
                    )
                )
    return stop_observations


def detect_stop_codons(
    path: Path,
    *,
    genetic_code: int | str | None = None,
) -> list[StopCodonObservation]:
    """Detect stop codons in a coding alignment under one chosen genetic code."""
    records = load_fasta_alignment(path)
    return _detect_stop_codons_in_records(records, genetic_code=genetic_code)


def _classify_sequence_coding_behavior_records(
    records: list[AlignmentRecord],
    *,
    genetic_code: int | str | None = None,
) -> list[SequenceCodingBehavior]:
    stop_observations = _detect_stop_codons_in_records(
        records, genetic_code=genetic_code
    )
    stop_counts_by_identifier: dict[str, list[StopCodonObservation]] = {}
    for stop in stop_observations:
        stop_counts_by_identifier.setdefault(stop.identifier, []).append(stop)
    invalid_observations = _detect_invalid_codons_in_records(records)
    invalid_counts_by_identifier: dict[str, list[InvalidCodonObservation]] = {}
    for invalid in invalid_observations:
        invalid_counts_by_identifier.setdefault(invalid.identifier, []).append(invalid)

    behaviors: list[SequenceCodingBehavior] = []
    for record in records:
        comparable_length = len(_coding_residues(record.sequence))
        stops = stop_counts_by_identifier.get(record.identifier, [])
        invalid_codons = invalid_counts_by_identifier.get(record.identifier, [])
        premature_stop_count = sum(1 for stop in stops if not stop.terminal)
        terminal_stop_count = sum(1 for stop in stops if stop.terminal)
        invalid_codon_count = len(invalid_codons)
        divisible_by_three = comparable_length % 3 == 0
        coding_like = (
            divisible_by_three
            and invalid_codon_count == 0
            and premature_stop_count == 0
        )
        if not comparable_length:
            note = "sequence contains no comparable coding residues after gaps and missing data are removed"
        elif coding_like and terminal_stop_count <= 1:
            note = "sequence is consistent with a coding reading frame"
        elif not divisible_by_three:
            note = (
                "sequence is not frame-consistent after removing gaps and missing data"
            )
        elif invalid_codon_count:
            note = "sequence contains ambiguous or invalid codons after normalization to coding triplets"
        elif premature_stop_count:
            note = "sequence contains one or more premature stop codons"
        else:
            note = "sequence shows mixed evidence for a coding interpretation"
        behaviors.append(
            SequenceCodingBehavior(
                identifier=record.identifier,
                coding_like=coding_like,
                comparable_length=comparable_length,
                divisible_by_three=divisible_by_three,
                invalid_codon_count=invalid_codon_count,
                premature_stop_count=premature_stop_count,
                terminal_stop_count=terminal_stop_count,
                note=note,
            )
        )
    return sorted(behaviors, key=lambda item: item.identifier)


def classify_sequence_coding_behavior(
    path: Path,
    *,
    genetic_code: int | str | None = None,
) -> list[SequenceCodingBehavior]:
    """Classify each sequence as coding-like or noncoding-like within a nucleotide dataset."""
    records = load_fasta_records(path)
    return _classify_sequence_coding_behavior_records(
        records,
        genetic_code=genetic_code,
    )


def _sequence_type_for_coding_preparation(
    path: Path,
    *,
    records: list[AlignmentRecord],
    sequence_type: AlignmentAlphabet | None,
) -> AlignmentAlphabet:
    if sequence_type is not None and sequence_type not in {"dna", "rna"}:
        raise InvalidAlignmentError(
            "codon-aware alignment requires dna or rna coding input"
        )
    detected = detect_fasta_sequence_type(path, records=records)
    effective = sequence_type if sequence_type is not None else detected.selected_type
    if effective not in {"dna", "rna"}:
        raise InvalidAlignmentError(
            f"codon-aware alignment requires nucleotide coding input: {detected.note}"
        )
    return effective


def prepare_coding_sequences_for_alignment(
    path: Path,
    *,
    sequence_type: AlignmentAlphabet | None = None,
    genetic_code: int | str | None = None,
) -> tuple[list[AlignmentRecord], CodingSequencePreparationReport]:
    """Filter raw coding sequences into a translation-ready set for codon-aware alignment."""
    records = load_fasta_records(path)
    genetic_code_id, genetic_code_name, _forward_table, _stop_codons = (
        _resolve_genetic_code_table(genetic_code)
    )
    effective_sequence_type = _sequence_type_for_coding_preparation(
        path,
        records=records,
        sequence_type=sequence_type,
    )
    frameshifts = {
        row.identifier: row for row in _detect_frameshift_like_records(records)
    }
    sequences_by_identifier = {record.identifier: record.sequence for record in records}

    accepted_records: list[AlignmentRecord] = []
    accepted_identifiers: list[str] = []
    excluded_sequences: list[CodingSequenceExclusion] = []
    invalid_codon_sequence_count = 0
    terminal_stop_sequence_count = 0

    for behavior in _classify_sequence_coding_behavior_records(
        records,
        genetic_code=genetic_code_id,
    ):
        normalized_sequence = _coding_residues(
            sequences_by_identifier[behavior.identifier]
        )
        if not normalized_sequence:
            excluded_sequences.append(
                CodingSequenceExclusion(
                    identifier=behavior.identifier,
                    comparable_length=0,
                    reason="empty-coding-sequence",
                    invalid_codon_count=0,
                    premature_stop_count=0,
                    terminal_stop_count=0,
                    trailing_bases=0,
                    note="sequence contains no comparable coding residues after gaps and missing data are removed",
                )
            )
            continue
        if not behavior.divisible_by_three:
            frameshift = frameshifts[behavior.identifier]
            excluded_sequences.append(
                CodingSequenceExclusion(
                    identifier=behavior.identifier,
                    comparable_length=behavior.comparable_length,
                    reason="frame-error",
                    invalid_codon_count=behavior.invalid_codon_count,
                    premature_stop_count=behavior.premature_stop_count,
                    terminal_stop_count=behavior.terminal_stop_count,
                    trailing_bases=frameshift.remainder,
                    note="sequence is not frame-consistent after gaps and missing data are removed",
                )
            )
            continue
        if behavior.invalid_codon_count:
            invalid_codon_sequence_count += 1
            excluded_sequences.append(
                CodingSequenceExclusion(
                    identifier=behavior.identifier,
                    comparable_length=behavior.comparable_length,
                    reason="invalid-codon",
                    invalid_codon_count=behavior.invalid_codon_count,
                    premature_stop_count=behavior.premature_stop_count,
                    terminal_stop_count=behavior.terminal_stop_count,
                    trailing_bases=0,
                    note="sequence contains ambiguous or invalid codons after normalization to coding triplets",
                )
            )
            continue
        if behavior.premature_stop_count:
            excluded_sequences.append(
                CodingSequenceExclusion(
                    identifier=behavior.identifier,
                    comparable_length=behavior.comparable_length,
                    reason="internal-stop-codon",
                    invalid_codon_count=behavior.invalid_codon_count,
                    premature_stop_count=behavior.premature_stop_count,
                    terminal_stop_count=behavior.terminal_stop_count,
                    trailing_bases=0,
                    note="sequence contains one or more premature stop codons",
                )
            )
            continue

        if behavior.terminal_stop_count:
            terminal_stop_sequence_count += 1
        normalized = normalized_sequence
        if effective_sequence_type == "rna":
            normalized = normalized.replace("T", "U")
        accepted_records.append(
            AlignmentRecord(identifier=behavior.identifier, sequence=normalized)
        )
        accepted_identifiers.append(behavior.identifier)

    if not accepted_records:
        raise InvalidAlignmentError(
            "codon-aware alignment excluded every sequence because of frame or stop-codon problems"
        )

    warnings: list[str] = []
    if excluded_sequences:
        warnings.append(
            "one or more coding sequences were excluded before codon-aware alignment"
        )
    if invalid_codon_sequence_count:
        warnings.append(
            "one or more coding sequences were excluded because they contained ambiguous or invalid codons"
        )
    if terminal_stop_sequence_count:
        warnings.append(
            "terminal stop codons were retained in accepted coding sequences"
        )
    return accepted_records, CodingSequencePreparationReport(
        source_path=path,
        sequence_type=effective_sequence_type,
        genetic_code_id=genetic_code_id,
        genetic_code_name=genetic_code_name,
        input_sequence_count=len(records),
        accepted_sequence_count=len(accepted_records),
        accepted_identifiers=accepted_identifiers,
        excluded_sequences=excluded_sequences,
        invalid_codon_sequence_count=invalid_codon_sequence_count,
        terminal_stop_sequence_count=terminal_stop_sequence_count,
        warnings=warnings,
    )


def inspect_coding_alignment(
    path: Path,
    *,
    genetic_code: int | str | None = None,
) -> CodingAlignmentDiagnostics:
    """Inspect one nucleotide alignment as a coding sequence dataset."""
    genetic_code_id, genetic_code_name, _forward_table, _stop_codons = (
        _resolve_genetic_code_table(genetic_code)
    )
    summary = summarise_fasta(path)
    if summary.inferred_alphabet not in {"dna", "rna"}:
        raise InvalidAlignmentError(
            f"coding diagnostics require a nucleotide alignment, got alphabet '{summary.inferred_alphabet}'"
        )
    records = load_fasta_alignment(path)
    coding_behaviors = _classify_sequence_coding_behavior_records(
        records,
        genetic_code=genetic_code_id,
    )
    invalid_codons = _detect_invalid_codons_in_records(records)
    partial_codon_sequences = [
        PartialCodonSequence(
            identifier=row.identifier,
            comparable_length=row.comparable_length,
            trailing_bases=row.remainder,
        )
        for row in _detect_frameshift_like_records(records)
    ]
    return CodingAlignmentDiagnostics(
        path=path,
        genetic_code_id=genetic_code_id,
        genetic_code_name=genetic_code_name,
        sequence_count=summary.sequence_count,
        alignment_length=summary.alignment_length,
        alignment_length_multiple_of_three=summary.alignment_length % 3 == 0,
        frameshift_like_sequences=_detect_frameshift_like_records(records),
        partial_codon_sequences=partial_codon_sequences,
        coding_behaviors=coding_behaviors,
        mixed_coding_signals=any(row.coding_like for row in coding_behaviors)
        and any(not row.coding_like for row in coding_behaviors),
        invalid_codons=invalid_codons,
        stop_codons=_detect_stop_codons_in_records(
            records,
            genetic_code=genetic_code_id,
        ),
    )


def _translate_codon(
    codon: str,
    *,
    genetic_code: int | str | None = None,
) -> str:
    _code_id, _code_name, forward_table, stop_codons = _resolve_genetic_code_table(
        genetic_code
    )
    normalized = codon.upper().replace("U", "T")
    if set(normalized) <= _GAP_CHARACTERS:
        return "-"
    if _invalid_codon_reason(normalized) is not None:
        return "X"
    if normalized in stop_codons:
        return "*"
    return forward_table.get(normalized, "X")


def _translation_status_for_codon(
    codon: str,
    amino_acid: str,
    *,
    codon_index: int,
    codon_count: int,
) -> str:
    normalized = codon.upper().replace("U", "T")
    if set(normalized) <= _GAP_CHARACTERS | _EXPLICIT_MISSING_CHARACTERS:
        return "missing-or-gap-codon"
    if _invalid_codon_reason(normalized) is not None:
        return "ambiguous-or-invalid-codon"
    if amino_acid == "*":
        if codon_index == codon_count:
            return "terminal-stop-codon"
        return "internal-stop-codon"
    return "translated"


def _build_translation_codon_observations(
    records: list[AlignmentRecord],
    *,
    translated_length: int,
    genetic_code_id: int,
) -> list[TranslationCodonObservation]:
    observations: list[TranslationCodonObservation] = []
    aligned_nucleotide_length = translated_length * 3
    for record in records:
        codons = _iter_codon_windows(record.sequence[:aligned_nucleotide_length])
        for codon_index, (nucleotide_start, codon) in enumerate(codons, start=1):
            amino_acid = _translate_codon(codon, genetic_code=genetic_code_id)
            observations.append(
                TranslationCodonObservation(
                    identifier=record.identifier,
                    codon_index=codon_index,
                    nucleotide_start=nucleotide_start,
                    codon=codon.upper(),
                    amino_acid=amino_acid,
                    translation_status=_translation_status_for_codon(
                        codon,
                        amino_acid,
                        codon_index=codon_index,
                        codon_count=translated_length,
                    ),
                )
            )
    return observations


def translate_coding_alignment(
    path: Path,
    *,
    genetic_code: int | str | None = None,
) -> tuple[list[AlignmentRecord], TranslationReport]:
    """Translate an aligned nucleotide coding sequence dataset to amino acids."""
    genetic_code_id, genetic_code_name, _forward_table, _stop_codons = (
        _resolve_genetic_code_table(genetic_code)
    )
    summary = summarise_fasta(path)
    if summary.inferred_alphabet not in {"dna", "rna"}:
        raise InvalidAlignmentError(
            f"coding translation requires a nucleotide alignment, got alphabet '{summary.inferred_alphabet}'"
        )

    records = load_fasta_alignment(path)
    dropped_trailing_nucleotide_count = summary.alignment_length % 3
    translated_alignment_length = (
        summary.alignment_length - dropped_trailing_nucleotide_count
    ) // 3
    warnings: list[str] = []
    if dropped_trailing_nucleotide_count:
        nucleotide_label = "nucleotide"
        if dropped_trailing_nucleotide_count != 1:
            nucleotide_label = "nucleotides"
        warnings.append(
            "sequence length not a multiple of 3: "
            f"{dropped_trailing_nucleotide_count} {nucleotide_label} dropped"
        )

    aligned_nucleotide_length = translated_alignment_length * 3
    translated_records = [
        AlignmentRecord(
            identifier=record.identifier,
            sequence="".join(
                _translate_codon(codon, genetic_code=genetic_code_id)
                for _, codon in _iter_codon_windows(
                    record.sequence[:aligned_nucleotide_length]
                )
            ),
        )
        for record in records
    ]
    codon_observations = _build_translation_codon_observations(
        records,
        translated_length=translated_alignment_length,
        genetic_code_id=genetic_code_id,
    )
    invalid_codon_count = sum(
        1
        for row in codon_observations
        if row.translation_status == "ambiguous-or-invalid-codon"
    )
    stop_codon_count = sum(1 for row in codon_observations if row.amino_acid == "*")
    internal_stop_sequence_count = len(
        {
            row.identifier
            for row in codon_observations
            if row.translation_status == "internal-stop-codon"
        }
    )
    terminal_stop_sequence_count = len(
        {
            row.identifier
            for row in codon_observations
            if row.translation_status == "terminal-stop-codon"
        }
    )
    return translated_records, TranslationReport(
        source_path=path,
        genetic_code_id=genetic_code_id,
        genetic_code_name=genetic_code_name,
        translated_sequence_count=len(translated_records),
        source_alignment_length=summary.alignment_length,
        translated_alignment_length=translated_alignment_length,
        dropped_trailing_nucleotide_count=dropped_trailing_nucleotide_count,
        invalid_codon_count=invalid_codon_count,
        stop_codon_count=stop_codon_count,
        internal_stop_sequence_count=internal_stop_sequence_count,
        terminal_stop_sequence_count=terminal_stop_sequence_count,
        trailing_partial_codon_sequence_count=(
            len(translated_records) if dropped_trailing_nucleotide_count else 0
        ),
        warnings=warnings,
        codon_observations=codon_observations,
        excluded_sequences=[],
    )


def write_translation_codon_validation_table(
    path: Path, report: TranslationReport
) -> Path:
    """Write a deterministic codon-validation ledger for one translation run."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "identifier\tcodon_index\tnucleotide_start\tcodon\tamino_acid\ttranslation_status"
    ]
    for row in report.codon_observations:
        lines.append(
            "\t".join(
                [
                    row.identifier,
                    str(row.codon_index),
                    str(row.nucleotide_start),
                    row.codon,
                    row.amino_acid,
                    row.translation_status,
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_translation_excluded_sequence_table(
    path: Path, report: TranslationReport
) -> Path:
    """Write a deterministic translation exclusion ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["identifier\treason\tnote"]
    for row in report.excluded_sequences:
        lines.append("\t".join([row.identifier, row.reason, row.note]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def translate_prepared_coding_sequences(
    records: list[AlignmentRecord],
    *,
    genetic_code: int | str | None = None,
) -> list[AlignmentRecord]:
    """Translate prepared coding sequences into an amino-acid guide alignment input."""
    genetic_code_id, _genetic_code_name, _forward_table, _stop_codons = (
        _resolve_genetic_code_table(genetic_code)
    )
    translated_records: list[AlignmentRecord] = []
    for record in records:
        if len(record.sequence) % 3 != 0:
            raise InvalidAlignmentError(
                "prepared coding sequences must remain divisible by three for translation"
            )
        translated_records.append(
            AlignmentRecord(
                identifier=record.identifier,
                sequence="".join(
                    "X" if amino_acid == "*" else amino_acid
                    for amino_acid in (
                        _translate_codon(codon, genetic_code=genetic_code_id)
                        for _, codon in _iter_codon_windows(record.sequence)
                    )
                ),
            )
        )
    return translated_records


def back_translate_aligned_coding_sequences(
    guide_alignment: list[AlignmentRecord],
    *,
    coding_records: list[AlignmentRecord],
) -> list[AlignmentRecord]:
    """Project an aligned amino-acid guide back onto nucleotide codon triplets."""
    coding_by_identifier = {record.identifier: record for record in coding_records}
    back_translated: list[AlignmentRecord] = []
    for guide_record in guide_alignment:
        try:
            coding_record = coding_by_identifier[guide_record.identifier]
        except KeyError as error:
            raise InvalidAlignmentError(
                f"codon back-translation is missing coding residues for {guide_record.identifier}"
            ) from error
        codons = [codon for _, codon in _iter_codon_windows(coding_record.sequence)]
        codon_index = 0
        aligned_codons: list[str] = []
        for residue in guide_record.sequence:
            if residue == "-":
                aligned_codons.append("---")
                continue
            if codon_index >= len(codons):
                raise InvalidAlignmentError(
                    f"aligned amino-acid guide consumed more codons than available for {guide_record.identifier}"
                )
            aligned_codons.append(codons[codon_index])
            codon_index += 1
        if codon_index != len(codons):
            raise InvalidAlignmentError(
                f"aligned amino-acid guide left unused codons for {guide_record.identifier}"
            )
        back_translated.append(
            AlignmentRecord(
                identifier=guide_record.identifier,
                sequence="".join(aligned_codons),
            )
        )
    return back_translated


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
    return _summarize_alignment_windows_from_records(
        summary,
        records,
        window_size=window_size,
        step_size=step_size,
    )


def _summarize_alignment_windows_from_records(
    summary: AlignmentSummary,
    records: list[AlignmentRecord],
    *,
    window_size: int,
    step_size: int,
) -> list[AlignmentWindowSummary]:
    windows: list[AlignmentWindowSummary] = []
    for start_index in range(0, summary.alignment_length, step_size):
        end_index = min(start_index + window_size, summary.alignment_length)
        if end_index <= start_index:
            continue
        # Flatten into per-position columns for metrics.
        columns = list(
            zip(
                *(record.sequence[start_index:end_index] for record in records),
                strict=True,
            )
        )
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
            missing_count += sum(
                1 for residue in column if _is_explicit_missing(residue)
            )
            ambiguity_count += sum(
                1
                for residue in column
                if _is_ambiguity_character(residue, alphabet=alphabet)
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
                state_counts = {
                    state: comparable.count(state) for state in set(comparable)
                }
                if len(state_counts) > 1:
                    variable_sites += 1
                disagreement_fractions.append(
                    1.0 - (max(state_counts.values()) / len(comparable))
                )
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
                disagreement_fraction=round(
                    sum(disagreement_fractions) / len(disagreement_fractions), 15
                ),
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
    return _detect_over_aligned_regions_from_windows(
        summarize_alignment_windows(path, window_size=window_size, step_size=step_size)
    )


def _detect_over_aligned_regions_from_windows(
    windows: list[AlignmentWindowSummary],
) -> list[AlignmentSuspiciousRegion]:
    regions: list[AlignmentSuspiciousRegion] = []
    for window in windows:
        uncertainty_fraction = (
            window.gap_fraction + window.missing_fraction + window.ambiguity_fraction
        )
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
    return _detect_under_aligned_regions_from_windows(
        summarize_alignment_windows(path, window_size=window_size, step_size=step_size)
    )


def _detect_under_aligned_regions_from_windows(
    windows: list[AlignmentWindowSummary],
) -> list[AlignmentSuspiciousRegion]:
    regions: list[AlignmentSuspiciousRegion] = []
    for window in windows:
        if window.variable_fraction >= 0.7 and window.disagreement_fraction >= 0.35:
            regions.append(
                AlignmentSuspiciousRegion(
                    start=window.start,
                    end=window.end,
                    kind="under_aligned",
                    score=round(
                        window.variable_fraction + window.disagreement_fraction, 15
                    ),
                    note="high local mismatch and disagreement suggest the region may require realignment or masking",
                )
            )
    return regions


def assess_alignment_low_information(
    path: Path,
    *,
    minimum_informative_sites: int = _LOW_INFORMATION_SITE_THRESHOLD,
    minimum_informative_fraction: float = _LOW_INFORMATION_FRACTION_THRESHOLD,
) -> AlignmentLowInformationReport:
    """Assess whether an alignment carries enough signal for defensible inference."""
    _validate_fraction_threshold(minimum_informative_fraction)
    if minimum_informative_sites < 0:
        raise ValueError(
            f"minimum_informative_sites must be non-negative, got {minimum_informative_sites}"
        )
    summary = summarise_fasta(path)
    return _assess_alignment_low_information_from_summary(
        summary,
        minimum_informative_sites=minimum_informative_sites,
        minimum_informative_fraction=minimum_informative_fraction,
    )


def _assess_alignment_low_information_from_summary(
    summary: AlignmentSummary,
    *,
    minimum_informative_sites: int,
    minimum_informative_fraction: float,
) -> AlignmentLowInformationReport:
    informative_fraction = (
        0.0
        if summary.alignment_length == 0
        else round(
            summary.parsimony_informative_site_count / summary.alignment_length,
            15,
        )
    )
    reasons: list[str] = []
    if (
        summary.parsimony_informative_site_count < minimum_informative_sites
        and summary.variable_site_count == 0
    ):
        reasons.append(
            "alignment has fewer parsimony-informative sites than the minimum threshold for defensible inference"
        )
    if (
        summary.variable_site_count == 0
        and informative_fraction < minimum_informative_fraction
    ):
        reasons.append(
            "alignment has a very low parsimony-informative-site fraction and may not support stable inference"
        )
    return AlignmentLowInformationReport(
        sequence_count=summary.sequence_count,
        alignment_length=summary.alignment_length,
        parsimony_informative_site_count=summary.parsimony_informative_site_count,
        parsimony_informative_fraction=informative_fraction,
        threshold_site_count=minimum_informative_sites,
        threshold_fraction=minimum_informative_fraction,
        low_information=bool(reasons),
        reasons=reasons,
    )


def build_duplicate_sequence_policy_report(
    path: Path,
    *,
    near_duplicate_threshold: float = 0.99,
) -> DuplicateSequencePolicyReport:
    """Build reviewer-facing policy recommendations for duplicate sequences."""
    exact_duplicates = detect_identical_duplicate_sequences(path)
    near_duplicates = detect_near_duplicate_sequences(
        path, identity_threshold=near_duplicate_threshold
    )
    actions: list[DuplicateSequencePolicyAction] = []
    warnings: list[str] = []
    if exact_duplicates:
        warnings.append(
            "exact duplicate sequences should be deduplicated or explicitly justified before inference"
        )
        for group in exact_duplicates:
            actions.append(
                DuplicateSequencePolicyAction(
                    action="collapse_exact_duplicates",
                    rationale=(
                        f"retain one representative such as {group.identifiers[0]} unless metadata show that "
                        "the duplicated labels represent distinct biological samples"
                    ),
                    affected_identifiers=group.identifiers,
                )
            )
    if near_duplicates:
        warnings.append(
            "near-duplicate sequences should be checked for replicate samples, contamination, or oversampling bias"
        )
        seen_pairs: set[tuple[str, ...]] = set()
        for pair in near_duplicates:
            identifiers = tuple(sorted((pair.left_identifier, pair.right_identifier)))
            if identifiers in seen_pairs:
                continue
            seen_pairs.add(identifiers)
            actions.append(
                DuplicateSequencePolicyAction(
                    action="review_near_duplicates",
                    rationale="inspect metadata and sampling provenance before keeping highly similar sequences together in inference",
                    affected_identifiers=list(identifiers),
                )
            )
    if not warnings:
        warnings.append("no duplicate-sequence policy actions are currently required")
    return DuplicateSequencePolicyReport(
        path=path,
        exact_duplicate_groups=exact_duplicates,
        near_duplicate_pairs=near_duplicates,
        policy_actions=actions,
        warnings=warnings,
    )


def build_ambiguous_alignment_column_report(
    path: Path,
    *,
    threshold: float = 0.5,
) -> AlignmentAmbiguousColumnReport:
    """Report columns dominated by ambiguity, explicit missingness, or gaps."""
    _validate_fraction_threshold(threshold)
    summary = summarise_fasta(path)
    return _build_ambiguous_alignment_column_report_from_summary(
        path,
        summary,
        threshold=threshold,
    )


def _build_ambiguous_alignment_column_report_from_summary(
    path: Path,
    summary: AlignmentSummary,
    *,
    threshold: float,
) -> AlignmentAmbiguousColumnReport:
    uncertainty_by_position = {
        row.position: row for row in summary.per_site_uncertainty
    }
    rows: list[AmbiguousAlignmentColumn] = []
    for position in range(1, summary.alignment_length + 1):
        uncertainty = uncertainty_by_position[position]
        ambiguity_burden = (
            uncertainty.ambiguity_fraction
            + uncertainty.missing_fraction
            + uncertainty.gap_fraction
        )
        if ambiguity_burden < threshold:
            continue
        rows.append(
            AmbiguousAlignmentColumn(
                position=position,
                ambiguity_fraction=uncertainty.ambiguity_fraction,
                missing_fraction=uncertainty.missing_fraction,
                gap_fraction=uncertainty.gap_fraction,
                comparable_fraction=round(max(0.0, 1.0 - ambiguity_burden), 15),
                note="column is dominated by ambiguity, missing data, or gaps and should be reviewed for masking",
            )
        )
    warnings = (
        [
            "alignment contains ambiguity-heavy columns that may be unsuitable for inference without masking"
        ]
        if rows
        else []
    )
    return AlignmentAmbiguousColumnReport(
        path=path,
        threshold=threshold,
        rows=rows,
        warnings=warnings,
    )


def build_sequence_quality_ranking(path: Path) -> SequenceQualityRankingReport:
    """Rank aligned sequences by transparent quality burdens."""
    summary = summarise_fasta(path)
    composition_outlier_ids = {row.identifier for row in summary.composition_outliers}
    exact_duplicate_ids = {
        identifier
        for group in summary.duplicate_sequence_groups
        for identifier in group.identifiers
    }
    near_duplicate_ids = {
        identifier
        for pair in summary.near_duplicate_pairs
        for identifier in (pair.left_identifier, pair.right_identifier)
    }
    uncertainty_by_id = {
        row.identifier: row for row in summary.per_sequence_uncertainty
    }
    ranked: list[tuple[float, str, SequenceQualityRankingRow]] = []
    for identifier in summary.ids:
        uncertainty = uncertainty_by_id[identifier]
        composition_outlier = identifier in composition_outlier_ids
        if identifier in exact_duplicate_ids:
            duplicate_status = "exact_duplicate"
        elif identifier in near_duplicate_ids:
            duplicate_status = "near_duplicate"
        else:
            duplicate_status = "unique"
        penalty = (
            uncertainty.missing_fraction * 40.0
            + uncertainty.gap_fraction * 25.0
            + uncertainty.ambiguity_fraction * 20.0
            + (10.0 if composition_outlier else 0.0)
            + (
                10.0
                if duplicate_status == "exact_duplicate"
                else 5.0
                if duplicate_status == "near_duplicate"
                else 0.0
            )
        )
        score = round(max(0.0, 100.0 - penalty), 3)
        note_parts: list[str] = []
        if uncertainty.missing_fraction > 0.0:
            note_parts.append("missing data")
        if uncertainty.gap_fraction > 0.0:
            note_parts.append("gaps")
        if uncertainty.ambiguity_fraction > 0.0:
            note_parts.append("ambiguity codes")
        if composition_outlier:
            note_parts.append("composition outlier")
        if duplicate_status != "unique":
            note_parts.append(duplicate_status.replace("_", " "))
        note = (
            "quality burdens: " + ", ".join(note_parts)
            if note_parts
            else "no major quality burdens detected"
        )
        ranked.append(
            (
                score,
                identifier,
                SequenceQualityRankingRow(
                    identifier=identifier,
                    rank=0,
                    score=score,
                    missing_fraction=uncertainty.missing_fraction,
                    gap_fraction=uncertainty.gap_fraction,
                    ambiguity_fraction=uncertainty.ambiguity_fraction,
                    composition_outlier=composition_outlier,
                    duplicate_status=duplicate_status,
                    note=note,
                ),
            )
        )
    ranked.sort(key=lambda item: (item[0], item[1]))
    rows = [
        SequenceQualityRankingRow(
            identifier=row.identifier,
            rank=index,
            score=row.score,
            missing_fraction=row.missing_fraction,
            gap_fraction=row.gap_fraction,
            ambiguity_fraction=row.ambiguity_fraction,
            composition_outlier=row.composition_outlier,
            duplicate_status=row.duplicate_status,
            note=row.note,
        )
        for index, (_, _, row) in enumerate(ranked, start=1)
    ]
    warnings = (
        ["lower-ranked sequences should be reviewed before publication or inference"]
        if rows and any(row.score < 85.0 for row in rows)
        else []
    )
    if not summary.near_duplicate_scan_performed:
        warnings.append(
            "near-duplicate sequence ranking was skipped because the alignment exceeds the governed pairwise review threshold"
        )
    return SequenceQualityRankingReport(path=path, rows=rows, warnings=warnings)


def summarize_alignment_readiness(path: Path) -> AlignmentReadinessReport:
    """Classify whether an input alignment is ready for key downstream analysis families."""
    sequence_kind = classify_alignment_sequences(path)
    records = load_fasta_records(path)
    inferred_alphabet = infer_alignment_alphabet(records)
    sequence_count = len(records)
    alignment_length = (
        sequence_kind.max_sequence_length
        if sequence_kind.state != "raw_sequence_fasta"
        else None
    )
    methods: list[AlignmentMethodReadiness] = []
    warnings: list[str] = []
    length_outliers = detect_sequence_length_outliers(path)

    if sequence_kind.state == "raw_sequence_fasta":
        common_blocker = ["input sequences are not yet aligned"]
        methods.extend(
            [
                AlignmentMethodReadiness(
                    analysis="distance",
                    ready=False,
                    blockers=common_blocker,
                    warnings=[],
                ),
                AlignmentMethodReadiness(
                    analysis="maximum_likelihood",
                    ready=False,
                    blockers=common_blocker,
                    warnings=[],
                ),
                AlignmentMethodReadiness(
                    analysis="bayesian",
                    ready=False,
                    blockers=common_blocker,
                    warnings=[],
                ),
                AlignmentMethodReadiness(
                    analysis="coding", ready=False, blockers=common_blocker, warnings=[]
                ),
                AlignmentMethodReadiness(
                    analysis="protein",
                    ready=inferred_alphabet == "protein",
                    blockers=[] if inferred_alphabet == "protein" else common_blocker,
                    warnings=[],
                ),
            ]
        )
        if length_outliers:
            warnings.append(
                "raw sequences include substantial length outliers before alignment"
            )
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
    low_information = assess_alignment_low_information(path)
    over_aligned = detect_over_aligned_regions(path)
    under_aligned = detect_under_aligned_regions(path)
    coding = (
        inspect_coding_alignment(path) if inferred_alphabet in {"dna", "rna"} else None
    )

    def _method(
        analysis: str, ready: bool, blockers: list[str], extra_warnings: list[str]
    ) -> AlignmentMethodReadiness:
        return AlignmentMethodReadiness(
            analysis=analysis, ready=ready, blockers=blockers, warnings=extra_warnings
        )

    generic_warnings: list[str] = []
    if sequence_kind.state == "ambiguous_equal_length_fasta":
        generic_warnings.append(
            "equal-length ungapped FASTA may be aligned, but prior alignment cannot be proved from sequence shape alone"
        )
    if over_aligned:
        generic_warnings.append("one or more windows look suspiciously over-aligned")
    if under_aligned:
        generic_warnings.append("one or more windows look suspiciously under-aligned")
    if quality.composition_outliers:
        generic_warnings.append("composition outliers may bias downstream inference")
    if quality.sequence_length_outliers:
        generic_warnings.append(
            "sequence length outliers suggest problematic trimming, concatenation, or mixed loci"
        )
    if coding is not None and coding.mixed_coding_signals:
        generic_warnings.append(
            "alignment mixes coding-like and noncoding-like sequence behavior"
        )
    if low_information.low_information:
        generic_warnings.extend(low_information.reasons)
    warnings.extend(generic_warnings)

    base_alignment_blockers = (
        ["alignment contains invalid characters for the inferred alphabet"]
        if quality.invalid_characters
        else []
    )
    distance_blockers = list(base_alignment_blockers)
    if summary.variable_site_count == 0:
        distance_blockers.append("alignment has no variable sites")
    if low_information.low_information:
        distance_blockers.append(
            "alignment has too few parsimony-informative sites for defensible inference"
        )
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
    if low_information.low_information:
        model_blockers.append(
            "alignment has too few parsimony-informative sites for defensible inference"
        )
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
            coding_blockers.append(
                "one or more sequences contain partial codons after gaps and missing data are removed"
            )
        if any(not stop.terminal for stop in coding.stop_codons):
            coding_blockers.append(
                "one or more sequences contain premature stop codons"
            )
        if any(stop.terminal for stop in coding.stop_codons):
            coding_warnings.append(
                "terminal stop codons were detected and should be verified against coding conventions"
            )
    methods.append(
        _method(
            "coding",
            ready=not coding_blockers,
            blockers=coding_blockers,
            extra_warnings=coding_warnings,
        )
    )

    protein_blockers = (
        []
        if inferred_alphabet == "protein"
        else ["protein analysis requires an amino-acid alignment"]
    )
    methods.append(
        _method(
            "protein",
            ready=not protein_blockers,
            blockers=protein_blockers,
            extra_warnings=generic_warnings,
        )
    )

    return AlignmentReadinessReport(
        path=path,
        sequence_kind=sequence_kind,
        inferred_alphabet=inferred_alphabet,
        sequence_count=summary.sequence_count,
        alignment_length=summary.alignment_length,
        methods=methods,
        warnings=warnings,
    )


def _alignment_quality_components(summary: AlignmentSummary) -> dict[str, float]:
    informative_density = (
        0.0
        if summary.alignment_length == 0
        else summary.parsimony_informative_site_count / summary.alignment_length
    )
    return {
        "missingness": round(
            max(0.0, 1.0 - min(summary.missing_data_fraction, 1.0)), 15
        ),
        "gap_burden": round(max(0.0, 1.0 - min(summary.gap_fraction, 1.0)), 15),
        "composition_outliers": round(
            max(
                0.0,
                1.0
                - min(
                    len(summary.composition_outliers) / max(summary.sequence_count, 1),
                    1.0,
                ),
            ),
            15,
        ),
        "duplicates": round(
            max(
                0.0,
                1.0
                - min(
                    len(summary.duplicate_sequence_groups)
                    / max(summary.sequence_count, 1),
                    1.0,
                ),
            ),
            15,
        ),
        "informative_density": round(min(informative_density / 0.2, 1.0), 15),
    }


def _alignment_quality_score(components: dict[str, float]) -> float:
    return round(sum(components.values()) / max(len(components), 1) * 100.0, 3)


def _summarize_missing_data_concentration(
    summary: AlignmentSummary,
    *,
    threshold: float = 0.5,
) -> AlignmentMissingDataConcentration:
    concentrated_columns = [
        row.position
        for row in summary.per_site_missingness
        if row.missing_fraction >= threshold
    ]
    longest_run = 0
    longest_start: int | None = None
    longest_end: int | None = None
    current_run = 0
    current_start: int | None = None
    previous_position: int | None = None
    for position in concentrated_columns:
        if previous_position is None or position != previous_position + 1:
            current_run = 1
            current_start = position
        else:
            current_run += 1
        if current_run > longest_run:
            longest_run = current_run
            longest_start = current_start
            longest_end = position
        previous_position = position
    maximum_missing_fraction = max(
        (row.missing_fraction for row in summary.per_site_missingness),
        default=0.0,
    )
    maximum_missing_positions = [
        row.position
        for row in summary.per_site_missingness
        if row.missing_fraction == maximum_missing_fraction
    ]
    return AlignmentMissingDataConcentration(
        threshold=threshold,
        concentrated_column_count=len(concentrated_columns),
        concentrated_column_fraction=(
            0.0
            if summary.alignment_length == 0
            else len(concentrated_columns) / summary.alignment_length
        ),
        longest_concentrated_run=longest_run,
        longest_concentrated_run_start=longest_start,
        longest_concentrated_run_end=longest_end,
        maximum_missing_fraction=maximum_missing_fraction,
        maximum_missing_positions=maximum_missing_positions,
    )


def _alignment_suspicion_reasons(
    *,
    low_information: AlignmentLowInformationReport,
    missing_data_concentration: AlignmentMissingDataConcentration,
    ambiguous_column_count: int,
    over_aligned_count: int,
    under_aligned_count: int,
    invalid_character_count: int,
) -> list[str]:
    reasons: list[str] = []
    if low_information.low_information:
        reasons.append("alignment has low information content for defensible inference")
    if missing_data_concentration.longest_concentrated_run >= 2:
        reasons.append("alignment concentrates missing data into adjacent columns")
    elif missing_data_concentration.concentrated_column_count > 0:
        reasons.append("alignment contains one or more highly missing columns")
    if ambiguous_column_count > 0:
        reasons.append("alignment contains ambiguity-heavy columns")
    if over_aligned_count > 0:
        reasons.append("alignment contains suspiciously over-aligned windows")
    if under_aligned_count > 0:
        reasons.append("alignment contains suspiciously under-aligned windows")
    if invalid_character_count > 0:
        reasons.append(
            "alignment contains invalid characters for the inferred alphabet"
        )
    return reasons


def build_alignment_quality_report(path: Path) -> AlignmentQualityReport:
    """Generate a higher-level alignment quality report from composition and identity diagnostics."""
    records = load_fasta_alignment(path)
    summary = summarise_records_as_alignment_summary(path=path, records=records)
    inferred_alphabet = summary.inferred_alphabet
    invalid_characters = summary.invalid_characters
    composition_outliers = summary.composition_outliers
    sequence_length_outliers = _detect_sequence_length_outlier_rows(
        [(record.identifier, len(record.sequence)) for record in records]
    )
    duplicate_sequence_groups = summary.duplicate_sequence_groups
    near_duplicate_pairs = summary.near_duplicate_pairs
    low_information = _assess_alignment_low_information_from_summary(
        summary,
        minimum_informative_sites=_LOW_INFORMATION_SITE_THRESHOLD,
        minimum_informative_fraction=_LOW_INFORMATION_FRACTION_THRESHOLD,
    )
    ambiguous_columns = _build_ambiguous_alignment_column_report_from_summary(
        path,
        summary,
        threshold=0.5,
    )
    windows = _summarize_alignment_windows_from_records(
        summary,
        records,
        window_size=30,
        step_size=10,
    )
    over_aligned = _detect_over_aligned_regions_from_windows(windows)
    under_aligned = _detect_under_aligned_regions_from_windows(windows)
    missing_data_concentration = _summarize_missing_data_concentration(summary)
    suspicious_reasons = _alignment_suspicion_reasons(
        low_information=low_information,
        missing_data_concentration=missing_data_concentration,
        ambiguous_column_count=len(ambiguous_columns.rows),
        over_aligned_count=len(over_aligned),
        under_aligned_count=len(under_aligned),
        invalid_character_count=len(invalid_characters),
    )
    quality_components = _alignment_quality_components(summary)
    warnings: list[str] = []
    if invalid_characters:
        warnings.append(
            "alignment contains characters invalid for the inferred alphabet"
        )
    if composition_outliers:
        warnings.append("alignment contains composition outlier sequences")
    if sequence_length_outliers:
        warnings.append("alignment contains raw-sequence length outliers")
    if duplicate_sequence_groups:
        warnings.append("alignment contains identical duplicate sequences")
    if near_duplicate_pairs:
        warnings.append("alignment contains near-duplicate sequences")
    if not summary.near_duplicate_scan_performed:
        warnings.append(
            "near-duplicate sequence scan was skipped because the alignment exceeds the governed pairwise review threshold"
        )
    warnings.extend(reason for reason in suspicious_reasons if reason not in warnings)
    return AlignmentQualityReport(
        path=path,
        sequence_count=summary.sequence_count,
        alignment_length=summary.alignment_length,
        invariant_site_count=summary.constant_site_count,
        missing_data_fraction=summary.missing_data_fraction,
        gap_fraction=summary.gap_fraction,
        ambiguity_fraction=summary.ambiguity_fraction,
        variable_site_count=summary.variable_site_count,
        parsimony_informative_site_count=summary.parsimony_informative_site_count,
        per_sequence_uncertainty=summary.per_sequence_uncertainty,
        per_site_uncertainty=summary.per_site_uncertainty,
        inferred_alphabet=inferred_alphabet,
        invalid_characters=invalid_characters,
        composition_outliers=composition_outliers,
        sequence_length_outliers=sequence_length_outliers,
        duplicate_sequence_groups=duplicate_sequence_groups,
        near_duplicate_pairs=near_duplicate_pairs,
        missing_data_concentration=missing_data_concentration,
        suspicious_alignment=bool(suspicious_reasons),
        suspicious_reasons=suspicious_reasons,
        quality_score=_alignment_quality_score(quality_components),
        quality_components=quality_components,
        warnings=warnings,
        near_duplicate_scan_performed=summary.near_duplicate_scan_performed,
    )


def build_alignment_forensic_report(path: Path) -> AlignmentForensicReport:
    """Integrate alignment quality, readiness, coding, and suspicious-region diagnostics."""
    quality = build_alignment_quality_report(path)
    readiness = summarize_alignment_readiness(path)
    summary = summarise_fasta(path)
    low_information = assess_alignment_low_information(path)
    coding = (
        inspect_coding_alignment(path)
        if summary.inferred_alphabet in {"dna", "rna"}
        else None
    )
    duplicate_policy = build_duplicate_sequence_policy_report(path)
    ambiguous_columns = build_ambiguous_alignment_column_report(path)
    sequence_ranking = build_sequence_quality_ranking(path)
    over_aligned = detect_over_aligned_regions(path)
    under_aligned = detect_under_aligned_regions(path)
    method_by_name = {method.analysis: method for method in readiness.methods}
    warnings = list(
        dict.fromkeys(
            [
                *quality.warnings,
                *readiness.warnings,
                *(
                    ["alignment mixes coding-like and noncoding-like sequences"]
                    if coding is not None and coding.mixed_coding_signals
                    else []
                ),
            ]
        )
    )
    limitations = [
        "alignment forensics summarize data quality and readiness but do not replace model-based tree inference validation",
        "publication readiness still depends on explicit method reporting, biological context, and reviewer inspection",
    ]
    return AlignmentForensicReport(
        path=path,
        quality=quality,
        readiness=readiness,
        low_information=low_information,
        coding=coding,
        duplicate_policy=duplicate_policy,
        ambiguous_columns=ambiguous_columns,
        sequence_ranking=sequence_ranking,
        over_aligned_regions=over_aligned,
        under_aligned_regions=under_aligned,
        safe_for_distance_analysis=method_by_name["distance"].ready,
        safe_for_maximum_likelihood=method_by_name["maximum_likelihood"].ready,
        safe_for_bayesian_inference=method_by_name["bayesian"].ready,
        safe_for_coding_analysis=method_by_name["coding"].ready,
        safe_for_publication=(
            quality.quality_score >= 75.0
            and not quality.suspicious_alignment
            and not warnings
        ),
        warnings=warnings,
        limitations=limitations,
    )


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
                identifier=current_identifier, sequence="".join(current_sequence)
            )
        )

    if not records:
        raise InvalidAlignmentError(f"alignment contains no FASTA records: {path}")

    return records


def classify_alignment_sequences(path: Path) -> AlignmentSequenceKindReport:
    """Classify whether a FASTA input already behaves like an aligned dataset."""
    records = load_permissive_fasta_records(path)
    lengths = [len(record.sequence) for record in records]
    min_length = min(lengths)
    max_length = max(lengths)
    has_gaps = any(
        any(residue in _GAP_CHARACTERS for residue in record.sequence)
        for record in records
    )

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
    records = load_permissive_fasta_records(path)
    return _detect_sequence_length_outlier_rows(
        [(record.identifier, len(record.sequence)) for record in records]
    )


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


def summarize_fasta_input(
    path: Path,
    *,
    records: list[AlignmentRecord] | None = None,
    sequence_type: AlignmentAlphabet | None = None,
) -> FastaInputSummary:
    """Summarize a FASTA input without assuming aligned equal-length sequences."""
    if records is None:
        scan = _scan_raw_fasta(path)
        sequence_type_report = _build_fasta_sequence_type_report(
            path=path,
            record_count=scan.record_count,
            shared_compatible=scan.shared_compatible,
            thymine_record_count=scan.thymine_record_count,
            uracil_record_count=scan.uracil_record_count,
            protein_signal_record_count=scan.protein_signal_record_count,
            invalid_record_count=scan.invalid_record_count,
        )
        return _build_fasta_input_summary_from_scan(
            scan,
            sequence_type=sequence_type,
            sequence_type_report=sequence_type_report,
        )
    lengths = [len(record.sequence) for record in records]
    sequence_type_report = detect_fasta_sequence_type(path, records=list(records))
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
        path=path,
        sequence_count=len(records),
        unique_identifier_count=len({record.identifier for record in records}),
        empty_sequence_count=sum(1 for record in records if not record.sequence),
        min_sequence_length=min(lengths),
        max_sequence_length=max(lengths),
        median_sequence_length=float(median(lengths)),
        total_residue_count=sum(lengths),
        inferred_alphabet=inferred_alphabet,
    )


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
        supported = (
            _DNA_CHARACTERS
            | _RNA_CHARACTERS
            | _PROTEIN_CHARACTERS
            | _GAP_CHARACTERS
            | _EXPLICIT_MISSING_CHARACTERS
        )
        allowed = supported
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


def validate_fasta_input(
    path: Path,
    *,
    sequence_type: AlignmentAlphabet | None = None,
) -> FastaInputValidationReport:
    """Validate a FASTA input before alignment or engine execution."""
    scan = _scan_raw_fasta(path)
    sequence_type_report = _build_fasta_sequence_type_report(
        path=path,
        record_count=scan.record_count,
        shared_compatible=scan.shared_compatible,
        thymine_record_count=scan.thymine_record_count,
        uracil_record_count=scan.uracil_record_count,
        protein_signal_record_count=scan.protein_signal_record_count,
        invalid_record_count=scan.invalid_record_count,
    )
    summary = _build_fasta_input_summary_from_scan(
        scan,
        sequence_type=sequence_type,
        sequence_type_report=sequence_type_report,
    )
    duplicate_identifiers = [
        FastaDuplicateIdentifier(
            identifier=identifier,
            occurrences=len(record_indices),
            record_indices=record_indices,
        )
        for identifier, record_indices in sorted(scan.duplicate_positions.items())
        if len(record_indices) > 1
    ]
    if summary.inferred_alphabet == "dna":
        illegal_characters = scan.illegal_dna
    elif summary.inferred_alphabet == "rna":
        illegal_characters = scan.illegal_rna
    elif summary.inferred_alphabet == "protein":
        illegal_characters = scan.illegal_protein
    else:
        illegal_characters = scan.illegal_supported
    empty_sequences = scan.empty_sequences
    length_outliers = _detect_sequence_length_outlier_rows(scan.lengths_by_identifier)
    warnings: list[str] = []
    if duplicate_identifiers:
        warnings.append("input contains duplicate sequence identifiers")
    if illegal_characters:
        warnings.append("input contains unsupported sequence characters")
    if empty_sequences:
        warnings.append("input contains empty sequences")
    if length_outliers:
        warnings.append("input contains sequence length outliers")
    warnings.extend(sequence_type_report.warnings)
    if (
        sequence_type is not None
        and sequence_type_report.selected_type is not None
        and sequence_type_report.selected_type != sequence_type
    ):
        warnings.append(
            "declared sequence type overrides automatic raw sequence classification"
        )
    if summary.inferred_alphabet == "unknown":
        warnings.append("input sequence type could not be inferred confidently")
    return FastaInputValidationReport(
        path=path,
        summary=summary,
        sequence_type_report=sequence_type_report,
        duplicate_identifiers=duplicate_identifiers,
        illegal_characters=illegal_characters,
        empty_sequences=empty_sequences,
        length_outliers=length_outliers,
        warnings=list(dict.fromkeys(warnings)),
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


def repair_fasta_input(
    path: Path,
    *,
    sequence_type: AlignmentAlphabet | None = None,
    normalize_identifiers: bool = False,
    remove_invalid_records: bool = False,
) -> tuple[list[AlignmentRecord], FastaRepairReport]:
    """Repair a FASTA input through explicit identifier and record policies."""
    records = load_permissive_fasta_records(path)
    validation = validate_fasta_input(path, sequence_type=sequence_type)
    if (
        sequence_type is None
        and validation.sequence_type_report.detected_type == "mixed"
    ):
        raise InvalidAlignmentError(
            "FASTA repair requires an explicit sequence_type when raw records carry conflicting sequence-type signals"
        )
    illegal_record_indices = {row.record_index for row in validation.illegal_characters}
    empty_record_indices = {row.record_index for row in validation.empty_sequences}
    mismatched_record_indices: set[int] = set()
    if remove_invalid_records and sequence_type is not None:
        for record_index, record in enumerate(records, start=1):
            compatible_types = _compatible_raw_sequence_types(record)
            if compatible_types and sequence_type not in compatible_types:
                mismatched_record_indices.add(record_index)
    retained_records: list[AlignmentRecord] = []
    normalized_identifiers: list[FastaIdentifierRepair] = []
    removed_records: list[FastaRemovedRecord] = []
    seen_identifiers: Counter[str] = Counter()

    for record_index, record in enumerate(records, start=1):
        removal_reasons: list[str] = []
        if remove_invalid_records and record_index in empty_record_indices:
            removal_reasons.append("empty-sequence")
        if remove_invalid_records and record_index in illegal_record_indices:
            removal_reasons.append("illegal-characters")
        if remove_invalid_records and record_index in mismatched_record_indices:
            removal_reasons.append("sequence-type-mismatch")
        if removal_reasons:
            removed_records.append(
                FastaRemovedRecord(
                    identifier=record.identifier,
                    record_index=record_index,
                    reason="+".join(removal_reasons),
                )
            )
            continue

        repaired_identifier = (
            _normalize_fasta_identifier(record.identifier)
            if normalize_identifiers
            else record.identifier
        )
        seen_identifiers[repaired_identifier] += 1
        duplicate_occurrence = seen_identifiers[repaired_identifier]
        if normalize_identifiers and duplicate_occurrence > 1:
            repaired_identifier = f"{repaired_identifier}_{duplicate_occurrence}"
        if repaired_identifier != record.identifier:
            note = "normalized identifier"
            if duplicate_occurrence > 1:
                note = "normalized identifier and resolved duplicate collision"
            normalized_identifiers.append(
                FastaIdentifierRepair(
                    original_identifier=record.identifier,
                    repaired_identifier=repaired_identifier,
                    record_index=record_index,
                    note=note,
                )
            )
        retained_records.append(
            AlignmentRecord(
                identifier=repaired_identifier,
                sequence=record.sequence,
            )
        )

    if not retained_records:
        raise InvalidAlignmentError(
            f"FASTA repair removed every sequence record from {path}"
        )

    after_summary = summarize_fasta_input(
        path,
        records=retained_records,
        sequence_type=sequence_type,
    )
    warnings = list(validation.warnings)
    if remove_invalid_records and removed_records:
        warnings.append("repair removed invalid sequence records")
    if mismatched_record_indices:
        warnings.append(
            "repair removed records incompatible with the declared sequence type"
        )
    if normalize_identifiers and normalized_identifiers:
        warnings.append("repair normalized FASTA identifiers")
    return retained_records, FastaRepairReport(
        source_path=path,
        output_path=None,
        before=validation.summary,
        after=after_summary,
        normalized_identifiers=normalized_identifiers,
        removed_records=removed_records,
        warnings=list(dict.fromkeys(warnings)),
    )


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


def summarise_records_as_alignment_summary(
    *, path: Path, records: list[AlignmentRecord]
) -> AlignmentSummary:
    """Summarise already-loaded equal-length alignment records."""
    inferred_alphabet = infer_alignment_alphabet(records)
    ids = [record.identifier for record in records]
    lengths = [len(record.sequence) for record in records]
    total_sites = len(records) * lengths[0]
    gap_count = sum(
        sum(1 for residue in record.sequence if residue in _GAP_CHARACTERS)
        for record in records
    )
    missing_count = sum(
        sum(
            1
            for residue in record.sequence
            if _is_missing_like(residue, alphabet=inferred_alphabet)
        )
        for record in records
    )
    ambiguity_count = sum(
        sum(
            1
            for residue in record.sequence
            if _is_ambiguity_character(residue, alphabet=inferred_alphabet)
        )
        for record in records
    )
    per_sequence_missingness = [
        SequenceMissingness(
            identifier=record.identifier,
            missing_fraction=sum(
                1
                for residue in record.sequence
                if _is_missing_like(residue, alphabet=inferred_alphabet)
            )
            / lengths[0],
        )
        for record in records
    ]
    per_sequence_uncertainty = [
        SequenceUncertaintyProfile(
            identifier=record.identifier,
            gap_fraction=sum(
                1 for residue in record.sequence if residue in _GAP_CHARACTERS
            )
            / lengths[0],
            missing_fraction=sum(
                1
                for residue in record.sequence
                if _is_missing_like(residue, alphabet=inferred_alphabet)
            )
            / lengths[0],
            ambiguity_fraction=sum(
                1
                for residue in record.sequence
                if _is_ambiguity_character(residue, alphabet=inferred_alphabet)
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
    for position, column in enumerate(
        zip(*(record.sequence for record in records), strict=True), start=1
    ):
        missing_fraction = sum(
            1
            for residue in column
            if _is_missing_like(residue, alphabet=inferred_alphabet)
        ) / len(records)
        per_site_missingness.append(
            SiteMissingness(position=position, missing_fraction=missing_fraction)
        )
        per_site_uncertainty.append(
            SiteUncertaintyProfile(
                position=position,
                gap_fraction=sum(1 for residue in column if residue in _GAP_CHARACTERS)
                / len(records),
                missing_fraction=missing_fraction,
                ambiguity_fraction=sum(
                    1
                    for residue in column
                    if _is_ambiguity_character(residue, alphabet=inferred_alphabet)
                )
                / len(records),
            )
        )
        if all(residue in _GAP_CHARACTERS for residue in column):
            all_gap_columns.append(position)
        if all(
            _is_missing_like(residue, alphabet=inferred_alphabet) for residue in column
        ):
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
    nucleotide_composition = compute_nucleotide_composition(
        records, alphabet=inferred_alphabet
    )
    amino_acid_composition = compute_amino_acid_composition(
        records, alphabet=inferred_alphabet
    )
    per_sequence_gc_content = compute_per_sequence_gc_content(
        records, alphabet=inferred_alphabet
    )
    whole_alignment_gc_content = compute_whole_alignment_gc_content(
        records, alphabet=inferred_alphabet
    )
    composition_outliers = _detect_composition_outlier_sequences_records(records)
    duplicate_sequence_groups = _detect_identical_duplicate_sequences_records(records)
    near_duplicate_pairs, near_duplicate_scan_performed = (
        _collect_near_duplicate_sequences_for_summary(
            records,
            identity_threshold=0.95,
        )
    )

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
        near_duplicate_scan_performed=near_duplicate_scan_performed,
    )


def summarise_fasta(path: Path) -> AlignmentSummary:
    """Summarise a FASTA alignment without loading a heavy dependency."""
    records = load_fasta_alignment(path)
    return summarise_records_as_alignment_summary(path=path, records=records)


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
    return [
        row
        for row in summary.per_sequence_missingness
        if row.missing_fraction > threshold
    ]


def detect_sites_with_excessive_missing_data(
    path: Path,
    *,
    threshold: float,
) -> list[SiteMissingness]:
    """Return alignment columns whose missing-data fraction exceeds the given threshold."""
    _validate_fraction_threshold(threshold)
    summary = summarise_fasta(path)
    return [
        row for row in summary.per_site_missingness if row.missing_fraction > threshold
    ]
