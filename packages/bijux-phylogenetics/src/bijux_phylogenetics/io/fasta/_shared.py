# ruff: noqa: F401
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from statistics import median

from Bio.Data import CodonTable

from bijux_phylogenetics.core.alignment import (
    AlignmentAlphabet,
    AlignmentAmbiguousColumnReport,
    AlignmentBaseFrequencyReport,
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
    AlignmentSegregatingSiteReport,
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
    DnaBinAlignment,
    DnaBinSequence,
    DnaBinStateRow,
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
    TrimmedAlignmentColumn,
)
from bijux_phylogenetics.core.metadata import load_taxon_table
from bijux_phylogenetics.runtime.errors import (
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
_APE_DNA_STATE_SET = set(_APE_DNA_STATE_ORDER)
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

def _normalize_ape_nucleotide_state(residue: str) -> str | None:
    normalized = residue.lower().replace("u", "t")
    if normalized in _APE_DNA_STATE_ORDER:
        return normalized
    return None

def _normalize_dnabin_state(
    residue: str,
    *,
    normalize_uracil: bool,
) -> str | None:
    normalized = residue.lower()
    if normalized == "u" and normalize_uracil:
        normalized = "t"
    if normalized in _APE_DNA_STATE_SET:
        return normalized
    return None

def _records_from_dnabin_alignment(
    alignment: DnaBinAlignment,
    *,
    uppercase: bool,
) -> list[AlignmentRecord]:
    return [
        AlignmentRecord(
            identifier=record.identifier,
            sequence=record.sequence.upper() if uppercase else record.sequence,
        )
        for record in alignment.records
    ]

def _ape_nucleotide_state_counts(sequence: str) -> dict[str, int]:
    counts = dict.fromkeys(_APE_DNA_STATE_ORDER, 0)
    for residue in sequence:
        normalized = _normalize_ape_nucleotide_state(residue)
        if normalized is None:
            continue
        counts[normalized] += 1
    return counts

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

def _ape_segregating_sequences_from_dna_bin_alignment(
    alignment: DnaBinAlignment,
) -> tuple[list[str], list[str]]:
    records = _records_from_dnabin_alignment(alignment, uppercase=False)
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
    return original_sequences, effective_sequences

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

def _composition_for_comparison(summary: AlignmentSummary) -> dict[str, float]:
    if summary.nucleotide_composition:
        return summary.nucleotide_composition
    return summary.amino_acid_composition

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
