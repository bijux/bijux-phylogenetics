from __future__ import annotations

from pathlib import Path
from statistics import median

from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    AlignmentRecord,
    FastaInputSummary,
    FastaSequenceTypeReport,
    SequenceLengthOutlier,
)

from .character_policy import (
    _PROTEIN_EXCLUSIVE_CHARACTERS_UPPER,
    _compatible_raw_sequence_types,
    _observed_raw_sequence_characters,
    _ordered_sequence_types,
)
from .raw_scan import _RawFastaScan, _scan_raw_fasta


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
