from __future__ import annotations

from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet as AlignmentAlphabet,
    AlignmentRecord as AlignmentRecord,
)

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
from .alignment_io import (
    detect_invalid_alignment_characters as detect_invalid_alignment_characters,
    load_fasta_alignment as load_fasta_alignment,
    load_fasta_records as load_fasta_records,
    load_permissive_fasta_records as load_permissive_fasta_records,
    write_fasta_alignment as write_fasta_alignment,
    _detect_invalid_alignment_characters_records as _detect_invalid_alignment_characters_records,
    _normalize_fasta_identifier as _normalize_fasta_identifier,
)

__all__ = [
    "AlignmentAlphabet",
    "AlignmentRecord",
    "compute_amino_acid_composition",
    "compute_nucleotide_composition",
    "compute_per_sequence_gc_content",
    "compute_whole_alignment_gc_content",
    "detect_fasta_sequence_type",
    "detect_invalid_alignment_characters",
    "infer_alignment_alphabet",
    "load_fasta_alignment",
    "load_fasta_records",
    "load_permissive_fasta_records",
    "write_fasta_alignment",
]
