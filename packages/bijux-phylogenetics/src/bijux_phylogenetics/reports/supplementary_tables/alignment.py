from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta.cleaning import compare_alignment_summaries
from bijux_phylogenetics.io.fasta.quality import (
    assess_alignment_low_information,
    build_alignment_quality_report,
    build_sequence_quality_ranking,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.phylo.alignment import (
    AlignmentLowInformationReport,
    AlignmentQualityReport,
    AlignmentSummary,
    SequenceQualityRankingReport,
    SequenceQualityRankingRow,
    SequenceUncertaintyProfile,
)

from .columns import alignment_table_columns
from .models import (
    SupplementaryAlignmentDiagnosticsRow,
    SupplementaryAlignmentDiagnosticsTableResult,
)
from .shared import stringify_list, write_dict_rows


def _serialize_alignment_row(
    row: SupplementaryAlignmentDiagnosticsRow,
) -> dict[str, object]:
    return {
        "sequence_id": row.sequence_id,
        "original_sequence_present": row.original_sequence_present,
        "filtered_sequence_present": row.filtered_sequence_present,
        "filtering_status": row.filtering_status,
        "filtering_reason": row.filtering_reason or "",
        "original_missing_fraction": ""
        if row.original_missing_fraction is None
        else row.original_missing_fraction,
        "original_gap_fraction": ""
        if row.original_gap_fraction is None
        else row.original_gap_fraction,
        "original_ambiguity_fraction": ""
        if row.original_ambiguity_fraction is None
        else row.original_ambiguity_fraction,
        "original_quality_score": ""
        if row.original_quality_score is None
        else row.original_quality_score,
        "duplicate_status": row.duplicate_status or "",
        "composition_outlier": ""
        if row.composition_outlier is None
        else row.composition_outlier,
        "original_alignment_length": ""
        if row.original_alignment_length is None
        else row.original_alignment_length,
        "original_sequence_count": ""
        if row.original_sequence_count is None
        else row.original_sequence_count,
        "original_missing_data_fraction": ""
        if row.original_missing_data_fraction is None
        else row.original_missing_data_fraction,
        "original_gap_fraction_alignment": ""
        if row.original_gap_fraction_alignment is None
        else row.original_gap_fraction_alignment,
        "original_variable_site_count": ""
        if row.original_variable_site_count is None
        else row.original_variable_site_count,
        "original_parsimony_informative_site_count": ""
        if row.original_parsimony_informative_site_count is None
        else row.original_parsimony_informative_site_count,
        "original_suspicious_alignment": ""
        if row.original_suspicious_alignment is None
        else row.original_suspicious_alignment,
        "original_low_information": ""
        if row.original_low_information is None
        else row.original_low_information,
        "original_low_information_reasons": stringify_list(
            row.original_low_information_reasons
        ),
        "filtered_alignment_length": ""
        if row.filtered_alignment_length is None
        else row.filtered_alignment_length,
        "filtered_sequence_count": ""
        if row.filtered_sequence_count is None
        else row.filtered_sequence_count,
        "filtered_missing_data_fraction": ""
        if row.filtered_missing_data_fraction is None
        else row.filtered_missing_data_fraction,
        "filtered_gap_fraction_alignment": ""
        if row.filtered_gap_fraction_alignment is None
        else row.filtered_gap_fraction_alignment,
        "filtered_variable_site_count": ""
        if row.filtered_variable_site_count is None
        else row.filtered_variable_site_count,
        "filtered_parsimony_informative_site_count": ""
        if row.filtered_parsimony_informative_site_count is None
        else row.filtered_parsimony_informative_site_count,
        "filtered_low_information": ""
        if row.filtered_low_information is None
        else row.filtered_low_information,
        "filtered_low_information_reasons": stringify_list(
            row.filtered_low_information_reasons
        ),
    }


def _sequence_uncertainty_lookup(
    summary: AlignmentSummary,
) -> dict[str, SequenceUncertaintyProfile]:
    return {row.identifier: row for row in summary.per_sequence_uncertainty}


def _sequence_quality_lookup(
    ranking: SequenceQualityRankingReport,
) -> dict[str, SequenceQualityRankingRow]:
    return {row.identifier: row for row in ranking.rows}


def _filtering_status(
    *,
    sequence_id: str,
    original_ids: set[str],
    filtered_ids: set[str] | None,
) -> tuple[str, str | None]:
    if filtered_ids is None:
        return "not_requested", None
    if sequence_id in original_ids and sequence_id in filtered_ids:
        return "retained_after_filtering", None
    if sequence_id in original_ids and sequence_id not in filtered_ids:
        return "removed_during_filtering", "absent_from_filtered_alignment"
    return "only_in_filtered_alignment", "absent_from_original_alignment"


def _alignment_sequence_order(
    original_summary: AlignmentSummary,
    filtered_summary: AlignmentSummary | None,
) -> list[str]:
    ordered = list(original_summary.ids)
    if filtered_summary is None:
        return ordered
    filtered_only = [
        identifier
        for identifier in filtered_summary.ids
        if identifier not in set(ordered)
    ]
    return [*ordered, *filtered_only]


def _build_alignment_row(
    *,
    sequence_id: str,
    original_summary: AlignmentSummary,
    original_quality: AlignmentQualityReport,
    original_low_information: AlignmentLowInformationReport,
    original_quality_lookup: dict[str, SequenceQualityRankingRow],
    original_uncertainty_lookup: dict[str, SequenceUncertaintyProfile],
    filtered_summary: AlignmentSummary | None,
    filtered_low_information: AlignmentLowInformationReport | None,
) -> SupplementaryAlignmentDiagnosticsRow:
    original_ids = set(original_summary.ids)
    filtered_ids = None if filtered_summary is None else set(filtered_summary.ids)
    filtering_status, filtering_reason = _filtering_status(
        sequence_id=sequence_id,
        original_ids=original_ids,
        filtered_ids=filtered_ids,
    )
    original_uncertainty = original_uncertainty_lookup.get(sequence_id)
    original_ranking = original_quality_lookup.get(sequence_id)
    return SupplementaryAlignmentDiagnosticsRow(
        sequence_id=sequence_id,
        original_sequence_present=sequence_id in original_ids,
        filtered_sequence_present=False
        if filtered_ids is None
        else sequence_id in filtered_ids,
        filtering_status=filtering_status,
        filtering_reason=filtering_reason,
        original_missing_fraction=None
        if original_uncertainty is None
        else original_uncertainty.missing_fraction,
        original_gap_fraction=None
        if original_uncertainty is None
        else original_uncertainty.gap_fraction,
        original_ambiguity_fraction=None
        if original_uncertainty is None
        else original_uncertainty.ambiguity_fraction,
        original_quality_score=None
        if original_ranking is None
        else original_ranking.score,
        duplicate_status=None
        if original_ranking is None
        else original_ranking.duplicate_status,
        composition_outlier=None
        if original_ranking is None
        else original_ranking.composition_outlier,
        original_alignment_length=original_summary.alignment_length,
        original_sequence_count=original_summary.sequence_count,
        original_missing_data_fraction=original_summary.missing_data_fraction,
        original_gap_fraction_alignment=original_summary.gap_fraction,
        original_variable_site_count=original_summary.variable_site_count,
        original_parsimony_informative_site_count=original_summary.parsimony_informative_site_count,
        original_suspicious_alignment=original_quality.suspicious_alignment,
        original_low_information=original_low_information.low_information,
        original_low_information_reasons=original_low_information.reasons,
        filtered_alignment_length=None
        if filtered_summary is None
        else filtered_summary.alignment_length,
        filtered_sequence_count=None
        if filtered_summary is None
        else filtered_summary.sequence_count,
        filtered_missing_data_fraction=None
        if filtered_summary is None
        else filtered_summary.missing_data_fraction,
        filtered_gap_fraction_alignment=None
        if filtered_summary is None
        else filtered_summary.gap_fraction,
        filtered_variable_site_count=None
        if filtered_summary is None
        else filtered_summary.variable_site_count,
        filtered_parsimony_informative_site_count=None
        if filtered_summary is None
        else filtered_summary.parsimony_informative_site_count,
        filtered_low_information=None
        if filtered_low_information is None
        else filtered_low_information.low_information,
        filtered_low_information_reasons=[]
        if filtered_low_information is None
        else filtered_low_information.reasons,
    )


def _write_alignment_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryAlignmentDiagnosticsRow],
) -> Path:
    return write_dict_rows(
        path,
        columns=columns,
        rows=[_serialize_alignment_row(row) for row in rows],
    )


def write_supplementary_alignment_diagnostics_table(
    path: Path,
    *,
    alignment_path: Path,
    filtered_alignment_path: Path | None = None,
) -> SupplementaryAlignmentDiagnosticsTableResult:
    """Write one supplementary alignment diagnostics table with optional filtering outcomes."""
    original_summary = summarise_fasta(alignment_path)
    original_quality = build_alignment_quality_report(alignment_path)
    original_low_information = assess_alignment_low_information(alignment_path)
    original_ranking = build_sequence_quality_ranking(alignment_path)
    filtered_summary = (
        None
        if filtered_alignment_path is None
        else summarise_fasta(filtered_alignment_path)
    )
    filtered_low_information = (
        None
        if filtered_alignment_path is None
        else assess_alignment_low_information(filtered_alignment_path)
    )
    if filtered_summary is not None:
        compare_alignment_summaries(alignment_path, original_summary, filtered_summary)
    sequence_ids = _alignment_sequence_order(original_summary, filtered_summary)
    original_uncertainty_lookup = _sequence_uncertainty_lookup(original_summary)
    original_quality_lookup = _sequence_quality_lookup(original_ranking)
    rows = [
        _build_alignment_row(
            sequence_id=sequence_id,
            original_summary=original_summary,
            original_quality=original_quality,
            original_low_information=original_low_information,
            original_quality_lookup=original_quality_lookup,
            original_uncertainty_lookup=original_uncertainty_lookup,
            filtered_summary=filtered_summary,
            filtered_low_information=filtered_low_information,
        )
        for sequence_id in sequence_ids
    ]
    columns = alignment_table_columns()
    _write_alignment_rows(path, columns=columns, rows=rows)
    retained_sequence_count = sum(
        1 for row in rows if row.filtering_status == "retained_after_filtering"
    )
    removed_sequence_count = sum(
        1 for row in rows if row.filtering_status == "removed_during_filtering"
    )
    filtered_only_sequence_count = sum(
        1 for row in rows if row.filtering_status == "only_in_filtered_alignment"
    )
    return SupplementaryAlignmentDiagnosticsTableResult(
        output_path=path,
        row_count=len(rows),
        retained_sequence_count=retained_sequence_count,
        removed_sequence_count=removed_sequence_count,
        filtered_only_sequence_count=filtered_only_sequence_count,
        columns=columns,
        rows=rows,
    )
