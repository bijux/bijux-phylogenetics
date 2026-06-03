from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.alignment import AlignmentQualityReport

from ..core import _detect_sequence_length_outlier_rows, load_fasta_alignment
from ..records import summarise_records_as_alignment_summary
from .site_diagnostics import (
    LOW_INFORMATION_FRACTION_THRESHOLD,
    LOW_INFORMATION_SITE_THRESHOLD,
    alignment_quality_components,
    alignment_quality_score,
    alignment_suspicion_reasons,
    assess_alignment_low_information_from_summary,
    build_ambiguous_alignment_column_report_from_summary,
    summarize_missing_data_concentration,
)
from .window_diagnostics import (
    detect_over_aligned_regions_from_windows,
    detect_under_aligned_regions_from_windows,
    summarize_alignment_windows_from_records,
)


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
    low_information = assess_alignment_low_information_from_summary(
        summary,
        minimum_informative_sites=LOW_INFORMATION_SITE_THRESHOLD,
        minimum_informative_fraction=LOW_INFORMATION_FRACTION_THRESHOLD,
    )
    ambiguous_columns = build_ambiguous_alignment_column_report_from_summary(
        path,
        summary,
        threshold=0.5,
    )
    windows = summarize_alignment_windows_from_records(
        summary,
        records,
        window_size=30,
        step_size=10,
    )
    over_aligned = detect_over_aligned_regions_from_windows(windows)
    under_aligned = detect_under_aligned_regions_from_windows(windows)
    missing_data_concentration = summarize_missing_data_concentration(summary)
    suspicious_reasons = alignment_suspicion_reasons(
        low_information=low_information,
        missing_data_concentration=missing_data_concentration,
        ambiguous_column_count=len(ambiguous_columns.rows),
        over_aligned_count=len(over_aligned),
        under_aligned_count=len(under_aligned),
        invalid_character_count=len(invalid_characters),
    )
    quality_components = alignment_quality_components(summary)
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
        quality_score=alignment_quality_score(quality_components),
        quality_components=quality_components,
        warnings=warnings,
        near_duplicate_scan_performed=summary.near_duplicate_scan_performed,
    )
