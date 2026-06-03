from __future__ import annotations

from pathlib import Path

from .models import MetadataJoinReport, MetadataJoinRow
from .tables import load_taxon_table


def join_table_to_taxa(
    taxa: list[str],
    path: Path,
    *,
    taxon_column: str | None = None,
) -> MetadataJoinReport:
    """Join a taxon-keyed table onto an explicit ordered taxon list."""
    table = load_taxon_table(path, taxon_column=taxon_column)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    joined_rows = [
        MetadataJoinRow(
            taxon=taxon,
            matched=taxon in rows_by_taxon,
            values=dict(rows_by_taxon[taxon]) if taxon in rows_by_taxon else {},
        )
        for taxon in taxa
    ]
    missing_from_metadata = [row.taxon for row in joined_rows if not row.matched]
    extra_metadata_taxa = sorted(set(table.taxa) - set(taxa))
    return MetadataJoinReport(
        path=path,
        tree_taxa=len(taxa),
        metadata_rows=table.row_count,
        taxon_column=table.taxon_column,
        joined_rows=joined_rows,
        missing_from_metadata=missing_from_metadata,
        extra_metadata_taxa=extra_metadata_taxa,
    )


__all__ = ["join_table_to_taxa"]
