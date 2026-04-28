from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.errors import MetadataJoinError


@dataclass(slots=True)
class TaxonTable:
    """Normalized metadata-style table keyed by a required taxon column."""

    path: Path
    format: str
    row_count: int
    columns: list[str]
    taxon_column: str
    rows: list[dict[str, str]]
    taxa: list[str]

    @property
    def indexed_values(self) -> set[str]:
        return set(self.taxa)

    @property
    def index_column(self) -> str:
        return self.taxon_column


@dataclass(slots=True)
class MetadataColumnCompleteness:
    """Completeness summary for one metadata column."""

    name: str
    missing_count: int
    completeness_fraction: float


@dataclass(slots=True)
class MetadataInspectionReport:
    """Stable summary of a metadata table keyed by taxon."""

    path: Path
    format: str
    row_count: int
    column_count: int
    columns: list[str]
    taxon_column: str
    taxa: list[str]
    column_completeness: list[MetadataColumnCompleteness]


def _detect_delimiter(path: Path) -> tuple[str, str]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return ",", "csv"
    if suffix == ".tsv":
        return "\t", "tsv"

    header_line = path.read_text(encoding="utf-8").splitlines()[0] if path.exists() else ""
    if "\t" in header_line:
        return "\t", "tsv"
    return ",", "csv"


def _resolve_taxon_column(columns: list[str], requested: str | None) -> str:
    if requested is not None:
        if requested not in columns:
            raise MetadataJoinError(f"missing taxon column '{requested}'")
        return requested

    preferred_columns = ("taxon", "taxa", "tip", "sample", "name", columns[0])
    for column in preferred_columns:
        if column in columns:
            return column
    raise MetadataJoinError("unable to resolve a taxon column")


def load_taxon_table(path: Path, *, taxon_column: str | None = None) -> TaxonTable:
    """Load a CSV or TSV table keyed by a unique taxon column."""
    if not path.exists():
        raise FileNotFoundError(f"table file not found: {path}")

    delimiter, table_format = _detect_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise MetadataJoinError(f"table has no header row: {path}")
        columns = [column.strip() for column in reader.fieldnames]
        resolved_taxon_column = _resolve_taxon_column(columns, taxon_column)

        rows: list[dict[str, str]] = []
        taxa: list[str] = []
        seen_taxa: set[str] = set()

        for row_index, row in enumerate(reader, start=2):
            normalized_row = {column: str(row.get(column, "")).strip() for column in columns}
            taxon = normalized_row[resolved_taxon_column]
            if not taxon:
                raise MetadataJoinError(
                    f"row {row_index} in {path} has an empty '{resolved_taxon_column}' value"
                )
            if taxon in seen_taxa:
                raise MetadataJoinError(f"duplicate taxon '{taxon}' found in {path}")
            seen_taxa.add(taxon)
            taxa.append(taxon)
            rows.append(normalized_row)

    return TaxonTable(
        path=path,
        format=table_format,
        row_count=len(rows),
        columns=columns,
        taxon_column=resolved_taxon_column,
        rows=rows,
        taxa=sorted(taxa),
    )


def inspect_metadata_table(path: Path, *, taxon_column: str | None = None) -> MetadataInspectionReport:
    """Inspect a metadata table and expose the stable taxon-key contract."""
    table = load_taxon_table(path, taxon_column=taxon_column)
    row_count = max(table.row_count, 1)
    column_completeness = [
        MetadataColumnCompleteness(
            name=column,
            missing_count=sum(1 for row in table.rows if not row[column]),
            completeness_fraction=sum(1 for row in table.rows if row[column]) / row_count,
        )
        for column in table.columns
    ]
    return MetadataInspectionReport(
        path=table.path,
        format=table.format,
        row_count=table.row_count,
        column_count=len(table.columns),
        columns=table.columns,
        taxon_column=table.taxon_column,
        taxa=table.taxa,
        column_completeness=column_completeness,
    )
