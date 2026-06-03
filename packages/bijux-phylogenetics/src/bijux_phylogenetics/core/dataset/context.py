from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.bayesian.beast.models import (
    FossilCalibrationValidationReport,
    TipDatingValidationReport,
)
from bijux_phylogenetics.bayesian.beast.validation import (
    validate_fossil_calibration_table,
    validate_tip_dating_metadata,
)
from bijux_phylogenetics.datasets.study_inputs import TaxonTable, load_taxon_table
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.fasta.quality import build_alignment_forensic_report
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment import AlignmentForensicReport

_GEOGRAPHY_COLUMN_HINTS = (
    "geography",
    "region",
    "location",
    "country",
    "area",
    "locality",
)
_EXTERNAL_ID_COLUMN_HINTS = (
    "accession",
    "accession_id",
    "sample_id",
    "isolate_id",
    "specimen_id",
    "ncbi_taxon_id",
    "gbif_id",
    "taxonomy_id",
    "external_id",
)
_GROUP_COLUMN_EXCLUSIONS = {"taxon", "taxa", "tip", "sample", "name"}


@dataclass(slots=True)
class _DatasetContext:
    tree_path: Path
    metadata_path: Path
    traits_path: Path
    alignment_path: Path | None
    tip_dates_path: Path | None
    calibration_path: Path | None
    tree_taxa: list[str]
    metadata_table: TaxonTable
    traits_table: TaxonTable
    alignment_ids: list[str]
    tip_date_taxa: list[str]
    geography_columns: list[str]
    geography_taxa: set[str]
    calibration_taxa_to_targets: dict[str, list[str]]
    alignment_forensic: AlignmentForensicReport | None
    tip_dates_report: TipDatingValidationReport | None
    calibration_report: FossilCalibrationValidationReport | None


def _ordered_taxa(table: TaxonTable) -> list[str]:
    return [row[table.taxon_column] for row in table.rows]


def _table_rows_by_taxon(table: TaxonTable) -> dict[str, dict[str, str]]:
    return {row[table.taxon_column]: row for row in table.rows}


def _is_external_id_column(name: str) -> bool:
    normalized = name.strip().lower()
    return normalized.endswith("_id") or normalized in _EXTERNAL_ID_COLUMN_HINTS


def _non_empty_taxa_for_columns(table: TaxonTable, columns: list[str]) -> set[str]:
    if not columns:
        return set()
    return {
        row[table.taxon_column]
        for row in table.rows
        if any(row[column].strip() for column in columns)
    }


def _detect_geography_columns(table: TaxonTable) -> list[str]:
    hints = set(_GEOGRAPHY_COLUMN_HINTS)
    return [
        column
        for column in table.columns
        if column != table.taxon_column and column.strip().lower() in hints
    ]


def _distinct_nonempty_values(table: TaxonTable, column: str) -> set[str]:
    return {row[column].strip() for row in table.rows if row[column].strip()}


def _infer_group_columns(table: TaxonTable) -> list[str]:
    columns: list[str] = []
    for column in table.columns:
        normalized = column.strip().lower()
        if normalized == table.taxon_column.strip().lower():
            continue
        if normalized in _GROUP_COLUMN_EXCLUSIONS or _is_external_id_column(normalized):
            continue
        values = _distinct_nonempty_values(table, column)
        if 2 <= len(values) <= 12:
            columns.append(column)
    return columns


def _calibration_taxa_to_targets(
    calibrations: FossilCalibrationValidationReport | None,
) -> dict[str, list[str]]:
    if calibrations is None:
        return {}
    mapping: dict[str, list[str]] = {}
    for calibration in calibrations.calibrations:
        if not calibration.valid:
            continue
        target = calibration.calibration_id
        for taxon in calibration.taxa:
            mapping.setdefault(taxon, []).append(target)
    return {taxon: sorted(targets) for taxon, targets in mapping.items()}


def _collect_external_ids(
    taxon: str,
    metadata_table: TaxonTable,
    traits_table: TaxonTable,
) -> dict[str, str]:
    external_ids: dict[str, str] = {}
    for surface_name, table in (("metadata", metadata_table), ("traits", traits_table)):
        row = _table_rows_by_taxon(table).get(taxon)
        if row is None:
            continue
        for column in table.columns:
            if column == table.taxon_column or not _is_external_id_column(column):
                continue
            value = row[column].strip()
            if value:
                external_ids[f"{surface_name}.{column}"] = value
    return external_ids


def _load_dataset_context(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
    alignment_forensic: AlignmentForensicReport | None = None,
    tip_dates_report: TipDatingValidationReport | None = None,
    calibration_report: FossilCalibrationValidationReport | None = None,
) -> _DatasetContext:
    tree = load_tree(tree_path)
    metadata_table = load_taxon_table(metadata_path)
    traits_table = load_taxon_table(traits_path)
    alignment_ids = (
        [record.identifier for record in load_fasta_alignment(alignment_path)]
        if alignment_path is not None
        else []
    )
    tip_dates_table = (
        load_taxon_table(tip_dates_path) if tip_dates_path is not None else None
    )
    geography_columns = sorted(
        set(_detect_geography_columns(metadata_table))
        | set(_detect_geography_columns(traits_table))
    )
    geography_taxa = _non_empty_taxa_for_columns(
        metadata_table, _detect_geography_columns(metadata_table)
    )
    geography_taxa.update(
        _non_empty_taxa_for_columns(
            traits_table, _detect_geography_columns(traits_table)
        )
    )
    calibrations = calibration_report
    if calibration_path is not None and calibrations is None:
        calibrations = validate_fossil_calibration_table(tree_path, calibration_path)
    tip_dates = tip_dates_report
    if tip_dates_path is not None and tip_dates is None:
        tip_dates = validate_tip_dating_metadata(
            tree_path, tip_dates_path, alignment_path=alignment_path
        )
    alignment_report = alignment_forensic
    if alignment_path is not None and alignment_report is None:
        alignment_report = build_alignment_forensic_report(alignment_path)
    return _DatasetContext(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
        tree_taxa=list(tree.tip_names),
        metadata_table=metadata_table,
        traits_table=traits_table,
        alignment_ids=alignment_ids,
        tip_date_taxa=[] if tip_dates_table is None else _ordered_taxa(tip_dates_table),
        geography_columns=geography_columns,
        geography_taxa=geography_taxa,
        calibration_taxa_to_targets=_calibration_taxa_to_targets(calibrations),
        alignment_forensic=alignment_report,
        tip_dates_report=tip_dates,
        calibration_report=calibrations,
    )
