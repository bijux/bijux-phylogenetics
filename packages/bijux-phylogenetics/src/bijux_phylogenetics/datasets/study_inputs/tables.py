from __future__ import annotations

from collections.abc import Mapping, Sequence
import csv
from pathlib import Path

from bijux_phylogenetics.runtime.errors import MetadataJoinError

from .models import TableValue, TaxonTable


def _detect_delimiter(path: Path) -> tuple[str, str]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return ",", "csv"
    if suffix == ".tsv":
        return "\t", "tsv"

    header_line = (
        path.read_text(encoding="utf-8").splitlines()[0] if path.exists() else ""
    )
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
            normalized_row = {
                column: _normalize_table_cell(row.get(column, "")) for column in columns
            }
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


def _normalize_table_cell(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _stringify_table_value(value: TableValue) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    return str(value)


def write_taxon_rows(
    path: Path, *, columns: list[str], rows: Sequence[Mapping[str, TableValue]]
) -> Path:
    """Write taxon-keyed rows as CSV or TSV based on the output suffix."""
    delimiter, _ = _detect_delimiter(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    stringified_rows = [
        {column: _stringify_table_value(row.get(column)) for column in columns}
        for row in rows
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(stringified_rows)
    return path


__all__ = ["TableValue", "TaxonTable", "load_taxon_table", "write_taxon_rows"]
