from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylo.alignment import (
    CodingSequencePreparationReport,
    SequenceCompositionOutlier,
)

from ..models import CatarrhineDataQualityStressPanelWorkflowReport
from .shared import _format_number, _substantive_alignment_warnings


def _write_raw_sequence_findings_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    rows: list[dict[str, str]] = []
    for row in report.raw_sequence_input_validation.duplicate_identifiers:
        rows.append(
            {
                "issue_kind": "duplicate_identifier",
                "identifier": row.identifier,
                "detail": (
                    f"occurrences={row.occurrences};record_indices="
                    + ",".join(str(value) for value in row.record_indices)
                ),
                "action": "normalize_identifier_collision",
            }
        )
    for row in report.raw_sequence_input_validation.illegal_characters:
        rows.append(
            {
                "issue_kind": "illegal_character",
                "identifier": row.identifier,
                "detail": f"record_index={row.record_index};position={row.position};character={row.character}",
                "action": "remove_invalid_record",
            }
        )
    for row in report.raw_sequence_input_validation.empty_sequences:
        rows.append(
            {
                "issue_kind": "empty_sequence",
                "identifier": row.identifier,
                "detail": f"record_index={row.record_index}",
                "action": "remove_invalid_record",
            }
        )
    for row in report.raw_sequence_length_outliers:
        rows.append(
            {
                "issue_kind": "sequence_length_outlier",
                "identifier": row.identifier,
                "detail": (
                    f"raw_length={row.raw_length};median_length="
                    f"{_format_number(row.median_length)};note={row.note}"
                ),
                "action": "drop_length_outlier",
            }
        )
    return write_taxon_rows(
        path,
        columns=["issue_kind", "identifier", "detail", "action"],
        rows=rows,
    )


def _write_raw_sequence_repair_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    rows: list[dict[str, str]] = []
    for row in report.raw_sequence_input_repair.normalized_identifiers:
        rows.append(
            {
                "action_kind": "rename_identifier",
                "identifier": row.original_identifier,
                "result_identifier": row.repaired_identifier,
                "detail": f"record_index={row.record_index};note={row.note}",
            }
        )
    for row in report.raw_sequence_input_repair.removed_records:
        rows.append(
            {
                "action_kind": "remove_record",
                "identifier": row.identifier,
                "result_identifier": "",
                "detail": f"record_index={row.record_index};reason={row.reason}",
            }
        )
    for row in report.raw_sequence_length_outliers:
        rows.append(
            {
                "action_kind": "drop_length_outlier",
                "identifier": row.identifier,
                "result_identifier": "",
                "detail": f"raw_length={row.raw_length};note={row.note}",
            }
        )
    return write_taxon_rows(
        path,
        columns=["action_kind", "identifier", "result_identifier", "detail"],
        rows=rows,
    )


def _write_repaired_sequence_validation_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    substantive_warnings = _substantive_alignment_warnings(
        report.repaired_sequence_input_validation.warnings
    )
    return write_taxon_rows(
        path,
        columns=[
            "surface",
            "sequence_count",
            "duplicate_identifier_count",
            "illegal_character_count",
            "empty_sequence_count",
            "length_outlier_count",
            "warning_count",
            "detail",
        ],
        rows=[
            {
                "surface": "repaired_sequence_input",
                "sequence_count": str(
                    report.repaired_sequence_input_validation.summary.sequence_count
                ),
                "duplicate_identifier_count": str(
                    len(report.repaired_sequence_input_validation.duplicate_identifiers)
                ),
                "illegal_character_count": str(
                    len(report.repaired_sequence_input_validation.illegal_characters)
                ),
                "empty_sequence_count": str(
                    len(report.repaired_sequence_input_validation.empty_sequences)
                ),
                "length_outlier_count": str(
                    len(report.repaired_sequence_input_validation.length_outliers)
                ),
                "warning_count": str(len(substantive_warnings)),
                "detail": (
                    "repaired FASTA input retains only unique identifiers and "
                    "records without illegal characters, empty bodies, or retained "
                    "length outliers"
                ),
            }
        ],
    )


def _write_coding_sequence_exclusions_table(
    path: Path,
    report: CodingSequencePreparationReport,
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "identifier",
            "reason",
            "comparable_length",
            "invalid_codon_count",
            "premature_stop_count",
            "terminal_stop_count",
            "trailing_bases",
            "note",
        ],
        rows=[
            {
                "identifier": row.identifier,
                "reason": row.reason,
                "comparable_length": str(row.comparable_length),
                "invalid_codon_count": str(row.invalid_codon_count),
                "premature_stop_count": str(row.premature_stop_count),
                "terminal_stop_count": str(row.terminal_stop_count),
                "trailing_bases": str(row.trailing_bases),
                "note": row.note,
            }
            for row in report.excluded_sequences
        ],
    )


def _write_sequence_outliers_table(
    path: Path,
    rows: list[SequenceCompositionOutlier],
    *,
    dropped_taxa: set[str],
) -> Path:
    return write_taxon_rows(
        path,
        columns=["taxon", "deviation", "robust_z_score", "action"],
        rows=[
            {
                "taxon": row.identifier,
                "deviation": _format_number(row.deviation),
                "robust_z_score": _format_number(row.robust_z_score),
                "action": (
                    "drop_taxon_from_cleaned_alignment"
                    if row.identifier in dropped_taxa
                    else "flag_only"
                ),
            }
            for row in rows
        ],
    )
