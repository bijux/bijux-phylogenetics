from __future__ import annotations

from ._shared import (
    AlignmentAlphabet,
    AlignmentLinkageReport,
    AlignmentRecord,
    AlignmentSequenceKindReport,
    AlignmentSummary,
    AlignmentTaxonMismatchError,
    Counter,
    FastaDuplicateIdentifier,
    FastaIdentifierRepair,
    FastaInputSummary,
    FastaInputValidationReport,
    FastaRemovedRecord,
    FastaRepairReport,
    InvalidAlignmentError,
    Path,
    SequenceLengthOutlier,
    SequenceMissingness,
    SequenceUncertaintyProfile,
    SiteMissingness,
    SiteUncertaintyProfile,
    _GAP_CHARACTERS,
    _build_fasta_input_summary_from_scan,
    _build_fasta_sequence_type_report,
    _collect_near_duplicate_sequences_for_summary,
    _compatible_raw_sequence_types,
    _detect_composition_outlier_sequences_records,
    _detect_identical_duplicate_sequences_records,
    _detect_sequence_length_outlier_rows,
    _is_ambiguity_character,
    _is_explicit_missing,
    _is_missing_like,
    _normalize_fasta_identifier,
    _scan_raw_fasta,
    _validate_fraction_threshold,
    compute_amino_acid_composition,
    compute_nucleotide_composition,
    compute_per_sequence_gc_content,
    compute_whole_alignment_gc_content,
    detect_fasta_sequence_type,
    detect_invalid_alignment_characters,
    infer_alignment_alphabet,
    load_fasta_alignment,
    load_permissive_fasta_records,
    load_tree,
    median,
)

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
