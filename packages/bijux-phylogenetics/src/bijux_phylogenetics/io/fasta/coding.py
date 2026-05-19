from __future__ import annotations

from ._shared import (
    AlignmentAlphabet,
    AlignmentRecord,
    CodingAlignmentDiagnostics,
    CodingSequenceExclusion,
    CodingSequencePreparationReport,
    DnaBinAlignment,
    FrameshiftLikeSequence,
    InvalidAlignmentError,
    PartialCodonSequence,
    Path,
    SequenceCodingBehavior,
    StopCodonObservation,
    TranslationReport,
    _build_translation_codon_observations,
    _classify_sequence_coding_behavior_records,
    _coding_residues,
    _detect_frameshift_like_records,
    _detect_invalid_codons_in_records,
    _detect_stop_codons_in_records,
    _iter_codon_windows,
    _records_from_dnabin_alignment,
    _resolve_genetic_code_table,
    _sequence_type_for_coding_preparation,
    _translate_codon,
    load_fasta_alignment,
    load_fasta_records,
)

from .matrix import (
    load_dna_bin_alignment,
)

def detect_frameshift_like_sequences(path: Path) -> list[FrameshiftLikeSequence]:
    """Detect coding sequences whose comparable length is not divisible by three."""
    records = load_fasta_alignment(path)
    return _detect_frameshift_like_records(records)


def detect_stop_codons(
    path: Path,
    *,
    genetic_code: int | str | None = None,
) -> list[StopCodonObservation]:
    """Detect stop codons in a coding alignment under one chosen genetic code."""
    records = load_fasta_alignment(path)
    return _detect_stop_codons_in_records(records, genetic_code=genetic_code)


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
    matrix = load_dna_bin_alignment(path, normalize_uracil=True)
    return inspect_coding_alignment_from_dna_bin_alignment(
        matrix,
        genetic_code=genetic_code,
    )


def inspect_coding_alignment_from_dna_bin_alignment(
    alignment: DnaBinAlignment,
    *,
    genetic_code: int | str | None = None,
) -> CodingAlignmentDiagnostics:
    """Inspect one DNAbin-compatible nucleotide matrix as a coding sequence dataset."""
    if alignment.source_alphabet not in {"dna", "rna"}:
        raise InvalidAlignmentError(
            f"coding diagnostics require a nucleotide alignment, got alphabet '{alignment.source_alphabet}'"
        )
    records = _records_from_dnabin_alignment(alignment, uppercase=True)
    genetic_code_id, genetic_code_name, _forward_table, _stop_codons = (
        _resolve_genetic_code_table(genetic_code)
    )
    coding_behaviors = _classify_sequence_coding_behavior_records(
        records,
        genetic_code=genetic_code_id,
    )
    partial_codon_sequences = [
        PartialCodonSequence(
            identifier=row.identifier,
            comparable_length=row.comparable_length,
            trailing_bases=row.remainder,
        )
        for row in _detect_frameshift_like_records(records)
    ]
    return CodingAlignmentDiagnostics(
        path=alignment.path,
        genetic_code_id=genetic_code_id,
        genetic_code_name=genetic_code_name,
        sequence_count=alignment.sequence_count,
        alignment_length=alignment.alignment_length,
        alignment_length_multiple_of_three=alignment.alignment_length % 3 == 0,
        frameshift_like_sequences=_detect_frameshift_like_records(records),
        partial_codon_sequences=partial_codon_sequences,
        coding_behaviors=coding_behaviors,
        mixed_coding_signals=any(row.coding_like for row in coding_behaviors)
        and any(not row.coding_like for row in coding_behaviors),
        invalid_codons=_detect_invalid_codons_in_records(records),
        stop_codons=_detect_stop_codons_in_records(
            records,
            genetic_code=genetic_code_id,
        ),
    )


def translate_coding_alignment(
    path: Path,
    *,
    genetic_code: int | str | None = None,
) -> tuple[list[AlignmentRecord], TranslationReport]:
    """Translate an aligned nucleotide coding sequence dataset to amino acids."""
    matrix = load_dna_bin_alignment(path, normalize_uracil=True)
    return translate_coding_alignment_from_dna_bin_alignment(
        matrix,
        genetic_code=genetic_code,
    )


def translate_coding_alignment_from_dna_bin_alignment(
    alignment: DnaBinAlignment,
    *,
    genetic_code: int | str | None = None,
) -> tuple[list[AlignmentRecord], TranslationReport]:
    """Translate one DNAbin-compatible coding matrix to amino acids."""
    if alignment.source_alphabet not in {"dna", "rna"}:
        raise InvalidAlignmentError(
            f"coding translation requires a nucleotide alignment, got alphabet '{alignment.source_alphabet}'"
        )
    genetic_code_id, genetic_code_name, _forward_table, _stop_codons = (
        _resolve_genetic_code_table(genetic_code)
    )
    records = _records_from_dnabin_alignment(alignment, uppercase=True)
    dropped_trailing_nucleotide_count = alignment.alignment_length % 3
    translated_alignment_length = (
        alignment.alignment_length - dropped_trailing_nucleotide_count
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
        source_path=alignment.path,
        genetic_code_id=genetic_code_id,
        genetic_code_name=genetic_code_name,
        translated_sequence_count=len(translated_records),
        source_alignment_length=alignment.alignment_length,
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
