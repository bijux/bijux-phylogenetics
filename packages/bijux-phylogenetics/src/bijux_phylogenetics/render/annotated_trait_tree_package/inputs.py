from __future__ import annotations

from bijux_phylogenetics.datasets.study_inputs import TaxonTable
from bijux_phylogenetics.render.tree_svg import AnnotationStrip
from bijux_phylogenetics.runtime.errors import MetadataJoinError


def require_table(
    table: TaxonTable | None, *, path: object | None, surface: str
) -> TaxonTable:
    """Require a concrete table for one annotated trait tree annotation surface."""
    if table is None or path is None:
        raise MetadataJoinError(
            f"annotated trait tree package requires {surface} input table"
        )
    return table


def build_full_label_map(
    *,
    taxa: list[str],
    metadata_table: TaxonTable | None,
    label_column: str | None,
) -> dict[str, str]:
    """Resolve the full reviewer-facing label map for rendered tree tips."""
    labels = {taxon: taxon for taxon in taxa}
    if label_column is None:
        return labels
    table = require_table(
        metadata_table,
        path=None if metadata_table is None else metadata_table.path,
        surface="a metadata table for label rendering",
    )
    if label_column not in table.columns:
        raise MetadataJoinError(
            f"metadata table does not contain label column '{label_column}'"
        )
    for row in table.rows:
        taxon = row[table.taxon_column]
        if row[label_column]:
            labels[taxon] = row[label_column]
    return labels


def build_string_map(table: TaxonTable, column: str) -> dict[str, str]:
    """Resolve one categorical annotation surface from a taxon table."""
    if column not in table.columns:
        raise MetadataJoinError(f"table does not contain column '{column}'")
    return {row[table.taxon_column]: row[column] for row in table.rows if row[column]}


def build_numeric_map(table: TaxonTable, column: str) -> dict[str, float]:
    """Resolve one numeric annotation surface from a taxon table."""
    if column not in table.columns:
        raise MetadataJoinError(f"table does not contain column '{column}'")
    values: dict[str, float] = {}
    for row in table.rows:
        raw_value = row[column]
        if not raw_value:
            continue
        try:
            values[row[table.taxon_column]] = float(raw_value)
        except ValueError as error:
            raise MetadataJoinError(
                f"column '{column}' contains a non-numeric value for taxon '{row[table.taxon_column]}'"
            ) from error
    return values


def build_annotation_strips(
    table: TaxonTable, columns: list[str]
) -> list[AnnotationStrip]:
    """Resolve named annotation strips from one supporting taxon table."""
    missing_columns = [column for column in columns if column not in table.columns]
    if missing_columns:
        raise MetadataJoinError(
            f"table does not contain columns: {', '.join(missing_columns)}"
        )
    return [
        AnnotationStrip(
            name=column,
            values={
                row[table.taxon_column]: row[column]
                for row in table.rows
                if row[column]
            },
        )
        for column in columns
    ]
