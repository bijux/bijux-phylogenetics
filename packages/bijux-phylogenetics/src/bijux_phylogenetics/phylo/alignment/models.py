from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

AlignmentAlphabet = str
AlignmentState = str
RawSequenceType = str
SequenceTypeConfidence = str


@dataclass(frozen=True, slots=True)
class AlignmentRecord:
    """Single FASTA alignment record."""

    identifier: str
    sequence: str


@dataclass(frozen=True, slots=True)
class AlignmentSequenceKindReport:
    """Classification of whether a FASTA input already behaves like an alignment."""

    path: Path
    sequence_count: int
    min_sequence_length: int
    max_sequence_length: int
    has_gap_characters: bool
    state: AlignmentState
    note: str


@dataclass(frozen=True, slots=True)
class FastaInputSummary:
    """Summary of a raw or aligned FASTA input before downstream validation."""

    path: Path
    sequence_count: int
    unique_identifier_count: int
    empty_sequence_count: int
    min_sequence_length: int
    max_sequence_length: int
    median_sequence_length: float
    total_residue_count: int
    inferred_alphabet: AlignmentAlphabet


@dataclass(frozen=True, slots=True)
class FastaSequenceTypeReport:
    """Automatic raw-sequence type classification for one FASTA input."""

    path: Path
    record_count: int
    detected_type: RawSequenceType
    selected_type: AlignmentAlphabet | None
    compatible_types: list[AlignmentAlphabet]
    confidence: SequenceTypeConfidence
    thymine_record_count: int
    uracil_record_count: int
    protein_signal_record_count: int
    invalid_record_count: int
    note: str
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class FastaDuplicateIdentifier:
    """One identifier observed more than once in a FASTA input."""

    identifier: str
    occurrences: int
    record_indices: list[int]


@dataclass(frozen=True, slots=True)
class FastaIllegalCharacter:
    """One unsupported sequence character observed in a FASTA input."""

    identifier: str
    record_index: int
    position: int
    character: str


@dataclass(frozen=True, slots=True)
class FastaEmptySequence:
    """One FASTA record whose sequence body is empty."""

    identifier: str
    record_index: int


@dataclass(frozen=True, slots=True)
class FastaIdentifierRepair:
    """One identifier rewrite applied during FASTA repair."""

    original_identifier: str
    repaired_identifier: str
    record_index: int
    note: str


@dataclass(frozen=True, slots=True)
class FastaRemovedRecord:
    """One record removed during FASTA repair."""

    identifier: str
    record_index: int
    reason: str


@dataclass(slots=True)
class FastaInputValidationReport:
    """Validation report for one FASTA sequence input before alignment."""

    path: Path
    summary: FastaInputSummary
    sequence_type_report: FastaSequenceTypeReport
    duplicate_identifiers: list[FastaDuplicateIdentifier]
    illegal_characters: list[FastaIllegalCharacter]
    empty_sequences: list[FastaEmptySequence]
    length_outliers: list[SequenceLengthOutlier]
    warnings: list[str]


@dataclass(slots=True)
class FastaRepairReport:
    """Repair report describing how a FASTA input changed."""

    source_path: Path
    output_path: Path | None
    before: FastaInputSummary
    after: FastaInputSummary
    normalized_identifiers: list[FastaIdentifierRepair]
    removed_records: list[FastaRemovedRecord]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class SequenceMissingness:
    """Missing-data fraction for one alignment sequence."""

    identifier: str
    missing_fraction: float


@dataclass(frozen=True, slots=True)
class SequenceUncertaintyProfile:
    """Per-sequence split of gaps, explicit missing data, and ambiguity codes."""

    identifier: str
    gap_fraction: float
    missing_fraction: float
    ambiguity_fraction: float


@dataclass(frozen=True, slots=True)
class SiteMissingness:
    """Missing-data fraction for one alignment column."""

    position: int
    missing_fraction: float


@dataclass(frozen=True, slots=True)
class SiteUncertaintyProfile:
    """Per-site split of gaps, explicit missing data, and ambiguity codes."""

    position: int
    gap_fraction: float
    missing_fraction: float
    ambiguity_fraction: float


@dataclass(frozen=True, slots=True)
class InvalidAlignmentCharacter:
    """One sequence character invalid for a declared alignment alphabet."""

    identifier: str
    position: int
    character: str


@dataclass(frozen=True, slots=True)
class SequenceGCContent:
    """GC content summary for one sequence."""

    identifier: str
    gc_fraction: float | None


@dataclass(frozen=True, slots=True)
class DnaBinStateRow:
    """One normalized DNAbin-compatible nucleotide state in an aligned matrix."""

    identifier: str
    position: int
    state: str


@dataclass(frozen=True, slots=True)
class DnaBinSequence:
    """One DNAbin-compatible aligned nucleotide sequence."""

    identifier: str
    sequence: str
    states: tuple[str, ...]


@dataclass(slots=True)
class DnaBinAlignment:
    """One stable DNAbin-compatible nucleotide alignment matrix."""

    path: Path
    source_alphabet: AlignmentAlphabet
    sequence_count: int
    alignment_length: int
    state_order: list[str]
    uracil_normalized: bool
    records: list[DnaBinSequence]
    rows: list[DnaBinStateRow]


@dataclass(frozen=True, slots=True)
class NucleotideStateFrequencyRow:
    """One reviewer-facing nucleotide state frequency row."""

    scope: str
    identifier: str | None
    state: str
    count: int
    frequency: float


@dataclass(frozen=True, slots=True)
class SegregatingSiteRow:
    """One reviewer-facing segregating alignment column."""

    position: int
    original_states: str
    effective_states: str
    known_state_count: int
    ambiguity_state_count: int
    gap_count: int
    missing_count: int


@dataclass(frozen=True, slots=True)
class SequenceCompositionOutlier:
    """One sequence whose composition deviates strongly from the alignment baseline."""

    identifier: str
    deviation: float
    robust_z_score: float | None


@dataclass(frozen=True, slots=True)
class SequenceLengthOutlier:
    """One raw sequence whose length is unusually short or long for the dataset."""

    identifier: str
    raw_length: int
    median_length: float
    robust_z_score: float | None
    note: str


@dataclass(frozen=True, slots=True)
class DuplicateSequenceGroup:
    """Identifiers sharing the exact same aligned sequence."""

    identifiers: list[str]
    sequence: str


@dataclass(frozen=True, slots=True)
class NearDuplicateSequencePair:
    """Pair of sequences above a caller-provided identity threshold."""

    left_identifier: str
    right_identifier: str
    identity: float
    comparable_sites: int


@dataclass(frozen=True, slots=True)
class TrimmedAlignmentColumn:
    """One removed alignment column with its 1-based position and reason."""

    position: int
    reason: str


@dataclass(frozen=True, slots=True)
class RemovedAlignmentSequence:
    """One removed alignment sequence with the explicit removal reason."""

    identifier: str
    missing_fraction: float
    reason: str


@dataclass(slots=True)
class AlignmentTrimReport:
    """Explicit report describing how an alignment was trimmed."""

    path: Path
    original_sequence_count: int
    trimmed_sequence_count: int
    original_alignment_length: int
    trimmed_alignment_length: int
    removed_columns: list[TrimmedAlignmentColumn]
    removed_sequences: list[RemovedAlignmentSequence]


@dataclass(frozen=True, slots=True)
class PairwiseSequenceIdentity:
    """Pairwise identity summary between two aligned sequences."""

    left_identifier: str
    right_identifier: str
    identity: float | None
    comparable_sites: int


@dataclass(slots=True)
class SequenceIdentityMatrix:
    """Deterministic matrix-style export of pairwise sequence identities."""

    path: Path
    identifiers: list[str]
    pairs: list[PairwiseSequenceIdentity]


@dataclass(frozen=True, slots=True)
class FrameshiftLikeSequence:
    """Sequence whose ungapped coding length is not divisible by three."""

    identifier: str
    comparable_length: int
    remainder: int


@dataclass(frozen=True, slots=True)
class PartialCodonSequence:
    """Sequence ending with an incomplete coding triplet after removing gaps and missing data."""

    identifier: str
    comparable_length: int
    trailing_bases: int


@dataclass(frozen=True, slots=True)
class SequenceCodingBehavior:
    """Coding-like versus noncoding-like behavior for one nucleotide sequence."""

    identifier: str
    coding_like: bool
    comparable_length: int
    divisible_by_three: bool
    invalid_codon_count: int
    premature_stop_count: int
    terminal_stop_count: int
    note: str


@dataclass(frozen=True, slots=True)
class CodingSequenceExclusion:
    """One raw coding sequence excluded before codon-aware alignment."""

    identifier: str
    comparable_length: int
    reason: str
    invalid_codon_count: int
    premature_stop_count: int
    terminal_stop_count: int
    trailing_bases: int
    note: str


@dataclass(frozen=True, slots=True)
class InvalidCodonObservation:
    """Observed invalid or ambiguous codon within a coding sequence."""

    identifier: str
    codon_index: int
    nucleotide_start: int
    codon: str
    reason: str


@dataclass(frozen=True, slots=True)
class StopCodonObservation:
    """Observed stop codon within an aligned coding sequence."""

    identifier: str
    codon_index: int
    nucleotide_start: int
    codon: str
    terminal: bool


@dataclass(slots=True)
class CodingAlignmentDiagnostics:
    """Coding-sequence diagnostics for one aligned nucleotide dataset."""

    path: Path
    genetic_code_id: int
    genetic_code_name: str
    sequence_count: int
    alignment_length: int
    alignment_length_multiple_of_three: bool
    frameshift_like_sequences: list[FrameshiftLikeSequence]
    partial_codon_sequences: list[PartialCodonSequence]
    coding_behaviors: list[SequenceCodingBehavior]
    mixed_coding_signals: bool
    invalid_codons: list[InvalidCodonObservation]
    stop_codons: list[StopCodonObservation]


@dataclass(slots=True)
class TranslationCodonObservation:
    """One codon-level observation from aligned nucleotide translation."""

    identifier: str
    codon_index: int
    nucleotide_start: int
    codon: str
    amino_acid: str
    translation_status: str


@dataclass(slots=True)
class TranslationSequenceExclusion:
    """One sequence excluded from aligned translation output."""

    identifier: str
    reason: str
    note: str


@dataclass(slots=True)
class TranslationReport:
    """Explicit record of a coding-alignment translation run."""

    source_path: Path
    genetic_code_id: int
    genetic_code_name: str
    translated_sequence_count: int
    source_alignment_length: int
    translated_alignment_length: int
    dropped_trailing_nucleotide_count: int
    invalid_codon_count: int
    stop_codon_count: int
    internal_stop_sequence_count: int
    terminal_stop_sequence_count: int
    trailing_partial_codon_sequence_count: int
    warnings: list[str]
    codon_observations: list[TranslationCodonObservation]
    excluded_sequences: list[TranslationSequenceExclusion]


@dataclass(slots=True)
class CodingSequencePreparationReport:
    """Preparation summary for codon-aware alignment from raw coding sequences."""

    source_path: Path
    sequence_type: AlignmentAlphabet
    genetic_code_id: int
    genetic_code_name: str
    input_sequence_count: int
    accepted_sequence_count: int
    accepted_identifiers: list[str]
    excluded_sequences: list[CodingSequenceExclusion]
    invalid_codon_sequence_count: int
    terminal_stop_sequence_count: int
    warnings: list[str]


@dataclass(slots=True)
class AlignmentQualityReport:
    """Higher-level alignment quality report built from composition and identity diagnostics."""

    path: Path
    sequence_count: int
    alignment_length: int
    invariant_site_count: int
    missing_data_fraction: float
    gap_fraction: float
    ambiguity_fraction: float
    variable_site_count: int
    parsimony_informative_site_count: int
    per_sequence_uncertainty: list[SequenceUncertaintyProfile]
    per_site_uncertainty: list[SiteUncertaintyProfile]
    inferred_alphabet: AlignmentAlphabet
    invalid_characters: list[InvalidAlignmentCharacter]
    composition_outliers: list[SequenceCompositionOutlier]
    sequence_length_outliers: list[SequenceLengthOutlier]
    duplicate_sequence_groups: list[DuplicateSequenceGroup]
    near_duplicate_pairs: list[NearDuplicateSequencePair]
    missing_data_concentration: AlignmentMissingDataConcentration
    suspicious_alignment: bool
    suspicious_reasons: list[str]
    quality_score: float
    quality_components: dict[str, float]
    warnings: list[str]
    near_duplicate_scan_performed: bool = True


@dataclass(frozen=True, slots=True)
class AlignmentMissingDataConcentration:
    """Compact summary of whether missing data clusters into concentrated columns."""

    threshold: float
    concentrated_column_count: int
    concentrated_column_fraction: float
    longest_concentrated_run: int
    longest_concentrated_run_start: int | None
    longest_concentrated_run_end: int | None
    maximum_missing_fraction: float
    maximum_missing_positions: list[int]


@dataclass(frozen=True, slots=True)
class AlignmentWindowSummary:
    """Sliding-window summary across an aligned FASTA dataset."""

    start: int
    end: int
    site_count: int
    gap_fraction: float
    missing_fraction: float
    ambiguity_fraction: float
    variable_fraction: float
    disagreement_fraction: float
    comparable_fraction: float


@dataclass(frozen=True, slots=True)
class AlignmentSuspiciousRegion:
    """One alignment window flagged as suspicious for over- or under-alignment."""

    start: int
    end: int
    kind: str
    score: float
    note: str


@dataclass(frozen=True, slots=True)
class AlignmentMethodReadiness:
    """Readiness decision for one downstream analysis family."""

    analysis: str
    ready: bool
    blockers: list[str]
    warnings: list[str]


@dataclass(slots=True)
class AlignmentReadinessReport:
    """Integrated readiness report for multiple downstream analysis families."""

    path: Path
    sequence_kind: AlignmentSequenceKindReport
    inferred_alphabet: AlignmentAlphabet
    sequence_count: int
    alignment_length: int | None
    methods: list[AlignmentMethodReadiness]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class AlignmentFilterProfile:
    """Named alignment-cleaning policy with explicit thresholds."""

    name: str
    remove_all_gap_sites: bool
    remove_all_missing_sites: bool
    site_missingness_threshold: float | None
    sequence_missingness_threshold: float | None
    preserve_codon_structure: bool
    note: str


@dataclass(frozen=True, slots=True)
class AlignmentCompositionShift:
    """One composition component before-versus-after cleaning or comparison."""

    component: str
    before: float
    after: float
    delta: float


@dataclass(frozen=True, slots=True)
class AlignmentComparisonReport:
    """Explicit comparison between two alignment versions."""

    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_only_taxa: list[str]
    right_only_taxa: list[str]
    left_alignment_length: int
    right_alignment_length: int
    left_missing_data_fraction: float
    right_missing_data_fraction: float
    left_gap_fraction: float
    right_gap_fraction: float
    left_variable_site_count: int
    right_variable_site_count: int
    left_parsimony_informative_site_count: int
    right_parsimony_informative_site_count: int
    composition_shifts: list[AlignmentCompositionShift]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class AlignmentSignalWarning:
    """Warning about loss of phylogenetic signal after cleaning."""

    code: str
    message: str


@dataclass(frozen=True, slots=True)
class AlignmentGroupRetention:
    """Retention summary for one metadata or trait group after cleaning."""

    column: str
    value: str
    original_count: int
    retained_count: int
    removed_count: int
    removed_fraction: float


@dataclass(frozen=True, slots=True)
class AlignmentLowInformationReport:
    """Assessment of whether an alignment contains enough phylogenetic signal for inference."""

    sequence_count: int
    alignment_length: int
    parsimony_informative_site_count: int
    parsimony_informative_fraction: float
    threshold_site_count: int
    threshold_fraction: float
    low_information: bool
    reasons: list[str]


@dataclass(frozen=True, slots=True)
class DuplicateSequencePolicyAction:
    """One reviewer-facing recommendation for handling duplicate or near-duplicate sequences."""

    action: str
    rationale: str
    affected_identifiers: list[str]


@dataclass(slots=True)
class DuplicateSequencePolicyReport:
    """Policy report describing how duplicate sequences should be handled before inference."""

    path: Path
    exact_duplicate_groups: list[DuplicateSequenceGroup]
    near_duplicate_pairs: list[NearDuplicateSequencePair]
    policy_actions: list[DuplicateSequencePolicyAction]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class AmbiguousAlignmentColumn:
    """One alignment column dominated by ambiguity or missing-like states."""

    position: int
    ambiguity_fraction: float
    missing_fraction: float
    gap_fraction: float
    comparable_fraction: float
    note: str


@dataclass(slots=True)
class AlignmentAmbiguousColumnReport:
    """Report of ambiguity-heavy columns that deserve masking or reviewer inspection."""

    path: Path
    threshold: float
    rows: list[AmbiguousAlignmentColumn]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class SequenceQualityRankingRow:
    """One ranked sequence with explicit quality burden components."""

    identifier: str
    rank: int
    score: float
    missing_fraction: float
    gap_fraction: float
    ambiguity_fraction: float
    composition_outlier: bool
    duplicate_status: str
    note: str


@dataclass(slots=True)
class SequenceQualityRankingReport:
    """Reviewer-facing ranking of sequence quality across one alignment."""

    path: Path
    rows: list[SequenceQualityRankingRow]
    warnings: list[str]


@dataclass(slots=True)
class AlignmentCleaningReport:
    """Profile-driven cleaning report with comparison and bias diagnostics."""

    profile: AlignmentFilterProfile
    trim: AlignmentTrimReport
    comparison: AlignmentComparisonReport
    signal_warnings: list[AlignmentSignalWarning]
    group_retention: list[AlignmentGroupRetention]
    warnings: list[str]


@dataclass(slots=True)
class AlignmentForensicReport:
    """Reviewer-facing summary of whether an alignment is safe for downstream use."""

    path: Path
    quality: AlignmentQualityReport
    readiness: AlignmentReadinessReport
    low_information: AlignmentLowInformationReport
    coding: CodingAlignmentDiagnostics | None
    duplicate_policy: DuplicateSequencePolicyReport
    ambiguous_columns: AlignmentAmbiguousColumnReport
    sequence_ranking: SequenceQualityRankingReport
    over_aligned_regions: list[AlignmentSuspiciousRegion]
    under_aligned_regions: list[AlignmentSuspiciousRegion]
    safe_for_distance_analysis: bool
    safe_for_maximum_likelihood: bool
    safe_for_bayesian_inference: bool
    safe_for_coding_analysis: bool
    safe_for_publication: bool
    warnings: list[str]
    limitations: list[str]


@dataclass(slots=True)
class AlignmentSummary:
    """Summary of an alignment input."""

    path: Path
    sequence_count: int
    alignment_length: int
    min_sequence_length: int
    max_sequence_length: int
    ids: list[str]
    missing_data_fraction: float
    gap_fraction: float
    ambiguity_fraction: float
    per_sequence_missingness: list[SequenceMissingness]
    per_sequence_uncertainty: list[SequenceUncertaintyProfile]
    per_site_missingness: list[SiteMissingness]
    per_site_uncertainty: list[SiteUncertaintyProfile]
    all_gap_columns: list[int]
    all_missing_columns: list[int]
    constant_site_count: int
    variable_site_count: int
    parsimony_informative_site_count: int
    inferred_alphabet: AlignmentAlphabet
    invalid_characters: list[InvalidAlignmentCharacter]
    nucleotide_composition: dict[str, float]
    amino_acid_composition: dict[str, float]
    per_sequence_gc_content: list[SequenceGCContent]
    whole_alignment_gc_content: float | None
    composition_outliers: list[SequenceCompositionOutlier]
    duplicate_sequence_groups: list[DuplicateSequenceGroup]
    near_duplicate_pairs: list[NearDuplicateSequencePair]
    near_duplicate_scan_performed: bool = True


@dataclass(slots=True)
class AlignmentBaseFrequencyReport:
    """ape-style nucleotide state frequencies for one alignment and its sequences."""

    path: Path
    inferred_alphabet: AlignmentAlphabet
    sequence_count: int
    alignment_length: int
    ambiguity_policy: str
    gap_policy: str
    missing_data_policy: str
    state_order: list[str]
    alignment_rows: list[NucleotideStateFrequencyRow]
    per_sequence_rows: list[NucleotideStateFrequencyRow]
    composition_outliers: list[SequenceCompositionOutlier]
    warnings: list[str]


@dataclass(slots=True)
class AlignmentSegregatingSiteReport:
    """ape-style segregating-site report for one nucleotide alignment."""

    path: Path
    inferred_alphabet: AlignmentAlphabet
    sequence_count: int
    alignment_length: int
    ambiguity_policy: str
    gap_policy: str
    missing_data_policy: str
    trailing_gap_policy: str
    segregating_site_positions: list[int]
    rows: list[SegregatingSiteRow]
    warnings: list[str]


@dataclass(slots=True)
class AlignmentLinkageReport:
    """Summary of how an alignment links against a tree tip set."""

    tree_path: Path
    alignment_path: Path
    tree_taxa: int
    alignment_ids: int
    linked_taxa: int
    usable_taxa: list[str]
    missing_from_alignment: list[str]
    extra_alignment_ids: list[str]
