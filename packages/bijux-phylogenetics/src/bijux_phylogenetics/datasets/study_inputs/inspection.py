from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.runtime.errors import MetadataJoinError

from .models import (
    MetadataColumnCompleteness,
    MetadataInspectionReport,
    TaxonTableIndexAudit,
)
from .tables import (
    _detect_delimiter,
    _normalize_table_cell,
    _resolve_taxon_column,
    load_taxon_table,
)


def inspect_taxon_table_index(
    path: Path, *, taxon_column: str | None = None
) -> TaxonTableIndexAudit:
    """Inspect taxon-key integrity without rejecting duplicate or empty keys early."""
    if not path.exists():
        raise FileNotFoundError(f"table file not found: {path}")

    delimiter, table_format = _detect_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise MetadataJoinError(f"table has no header row: {path}")
        columns = [column.strip() for column in reader.fieldnames]
        resolved_taxon_column = _resolve_taxon_column(columns, taxon_column)
        observed_taxa: list[str] = []
        duplicate_taxa: set[str] = set()
        empty_taxon_rows: list[int] = []
        seen_taxa: set[str] = set()
        row_count = 0
        for row_index, row in enumerate(reader, start=2):
            row_count += 1
            taxon = _normalize_table_cell(row.get(resolved_taxon_column, ""))
            if not taxon:
                empty_taxon_rows.append(row_index)
                continue
            observed_taxa.append(taxon)
            if taxon in seen_taxa:
                duplicate_taxa.add(taxon)
                continue
            seen_taxa.add(taxon)
    return TaxonTableIndexAudit(
        path=path,
        format=table_format,
        row_count=row_count,
        taxon_column=resolved_taxon_column,
        taxa=sorted(set(observed_taxa)),
        duplicate_taxa=sorted(duplicate_taxa),
        empty_taxon_rows=empty_taxon_rows,
    )


def inspect_metadata_table(
    path: Path, *, taxon_column: str | None = None
) -> MetadataInspectionReport:
    """Inspect a metadata table and expose the stable taxon-key contract."""
    table = load_taxon_table(path, taxon_column=taxon_column)
    row_count = max(table.row_count, 1)
    column_completeness = [
        MetadataColumnCompleteness(
            name=column,
            missing_count=sum(1 for row in table.rows if not row[column]),
            completeness_fraction=sum(1 for row in table.rows if row[column])
            / row_count,
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


__all__ = ["inspect_metadata_table", "inspect_taxon_table_index"]
