from __future__ import annotations

from ._shared import (
    AlignmentBaseFrequencyReport,
    AlignmentSegregatingSiteReport,
    DnaBinAlignment,
    DnaBinSequence,
    DnaBinStateRow,
    InvalidAlignmentError,
    NucleotideStateFrequencyRow,
    Path,
    SegregatingSiteRow,
    _APE_DNA_STATE_ORDER,
    _APE_SEGREGATING_STATE_SETS,
    _ape_nucleotide_state_counts,
    _ape_segregating_sequences_from_dna_bin_alignment,
    _detect_composition_outlier_sequences_records,
    _is_ape_gap_state,
    _is_ape_known_base,
    _is_ape_missing_state,
    _is_ape_segregating_column,
    _normalize_dnabin_state,
    _records_from_dnabin_alignment,
    infer_alignment_alphabet,
    load_fasta_alignment,
    write_fasta_alignment,
)

def load_dna_bin_alignment(
    path: Path,
    *,
    normalize_uracil: bool = False,
) -> DnaBinAlignment:
    """Load one equal-length nucleotide FASTA input as a DNAbin-compatible matrix."""
    records = load_fasta_alignment(path)
    source_alphabet = infer_alignment_alphabet(records)
    if source_alphabet == "unknown":
        source_alphabet = "dna"

    invalid_rows: list[tuple[str, int, str]] = []
    dnabin_records: list[DnaBinSequence] = []
    rows: list[DnaBinStateRow] = []
    for record in records:
        normalized_states: list[str] = []
        for position, residue in enumerate(record.sequence, start=1):
            normalized = _normalize_dnabin_state(
                residue,
                normalize_uracil=normalize_uracil,
            )
            if normalized is None:
                invalid_rows.append((record.identifier, position, residue))
                continue
            normalized_states.append(normalized)
            rows.append(
                DnaBinStateRow(
                    identifier=record.identifier,
                    position=position,
                    state=normalized,
                )
            )
        if len(normalized_states) != len(record.sequence):
            continue
        dnabin_records.append(
            DnaBinSequence(
                identifier=record.identifier,
                sequence="".join(normalized_states),
                states=tuple(normalized_states),
            )
        )

    if invalid_rows:
        details = ", ".join(
            f"{identifier}:{position}={character}"
            for identifier, position, character in invalid_rows[:5]
        )
        if len(invalid_rows) > 5:
            details = f"{details}, ..."
        raise InvalidAlignmentError(
            "alignment contains states incompatible with dnabin-compatible nucleotide loading: "
            + details
        )
    if source_alphabet not in {"dna", "rna"}:
        raise InvalidAlignmentError(
            f"dnabin-compatible loading requires a dna or rna alignment, got alphabet '{source_alphabet}'"
        )

    return DnaBinAlignment(
        path=path,
        source_alphabet=source_alphabet,
        sequence_count=len(dnabin_records),
        alignment_length=0 if not dnabin_records else len(dnabin_records[0].sequence),
        state_order=list(_APE_DNA_STATE_ORDER),
        uracil_normalized=normalize_uracil and source_alphabet == "rna",
        records=dnabin_records,
        rows=rows,
    )


def write_dna_bin_alignment_fasta(path: Path, alignment: DnaBinAlignment) -> Path:
    """Write one DNAbin-compatible nucleotide matrix back to FASTA without state loss."""
    return write_fasta_alignment(
        path,
        _records_from_dnabin_alignment(alignment, uppercase=False),
    )


def compute_alignment_base_frequency_report_from_dna_bin_alignment(
    alignment: DnaBinAlignment,
) -> AlignmentBaseFrequencyReport:
    """Compute ape-style nucleotide state frequencies from one DNAbin-compatible matrix."""
    records = _records_from_dnabin_alignment(alignment, uppercase=False)
    alphabet = alignment.source_alphabet

    alignment_counts = dict.fromkeys(_APE_DNA_STATE_ORDER, 0)
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
    canonical_total = sum(alignment_counts[state] for state in _APE_DNA_STATE_ORDER[:4])
    if canonical_total == 0:
        warnings.append(
            "alignment contains no canonical A/C/G/T residues, so ape-style base frequencies reflect only ambiguity, gap, and missing states"
        )
    return AlignmentBaseFrequencyReport(
        path=alignment.path,
        inferred_alphabet=alphabet,
        sequence_count=alignment.sequence_count,
        alignment_length=alignment.alignment_length,
        ambiguity_policy="count ambiguity codes as literal states",
        gap_policy="count gap characters as literal states",
        missing_data_policy="count explicit missing characters as literal states",
        state_order=list(_APE_DNA_STATE_ORDER),
        alignment_rows=alignment_rows,
        per_sequence_rows=per_sequence_rows,
        composition_outliers=_detect_composition_outlier_sequences_records(records),
        warnings=warnings,
    )


def compute_alignment_base_frequency_report(path: Path) -> AlignmentBaseFrequencyReport:
    """Compute ape-style nucleotide state frequencies for one DNA or RNA alignment."""
    matrix = load_dna_bin_alignment(path, normalize_uracil=True)
    return compute_alignment_base_frequency_report_from_dna_bin_alignment(matrix)


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


def compute_alignment_segregating_site_report_from_dna_bin_alignment(
    alignment: DnaBinAlignment,
) -> AlignmentSegregatingSiteReport:
    """Compute ape-style segregating sites from one DNAbin-compatible matrix."""
    alphabet = alignment.source_alphabet
    original_sequences, effective_sequences = (
        _ape_segregating_sequences_from_dna_bin_alignment(alignment)
    )

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
                gap_count=sum(
                    1 for state in effective_states if _is_ape_gap_state(state)
                ),
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
        path=alignment.path,
        inferred_alphabet=alphabet,
        sequence_count=alignment.sequence_count,
        alignment_length=alignment.alignment_length,
        ambiguity_policy="ambiguity states segregate only when they are surely incompatible with another observed state",
        gap_policy="internal gap characters can create segregating sites against known or incompatible ambiguous states",
        missing_data_policy="explicit missing characters do not create segregating sites",
        trailing_gap_policy="leading and trailing gap runs are normalized to N before ape-style segregating-site detection",
        segregating_site_positions=segregating_site_positions,
        rows=rows,
        warnings=warnings,
    )


def compute_alignment_segregating_site_report(
    path: Path,
) -> AlignmentSegregatingSiteReport:
    """Compute ape-style segregating sites for one DNA or RNA alignment."""
    matrix = load_dna_bin_alignment(path, normalize_uracil=True)
    return compute_alignment_segregating_site_report_from_dna_bin_alignment(matrix)


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
