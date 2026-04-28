from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.metadata import TaxonTable, load_taxon_table


@dataclass(slots=True)
class TraitColumnSummary:
    """Deterministic schema summary for one trait column."""

    name: str
    kind: str
    missing_count: int
    distinct_value_count: int


@dataclass(slots=True)
class TraitValidationReport:
    """Stable summary of a validated trait table."""

    path: Path
    format: str
    row_count: int
    taxon_column: str
    trait_columns: list[TraitColumnSummary]


def load_tsv_summary(path: Path) -> TaxonTable:
    """Compatibility wrapper for legacy callers expecting a taxon-keyed table."""
    return load_taxon_table(path)


def _is_numeric(values: list[str]) -> bool:
    try:
        for value in values:
            float(value)
    except ValueError:
        return False
    return True


def _summarize_trait_column(table: TaxonTable, column: str) -> TraitColumnSummary:
    values = [row[column] for row in table.rows]
    observed_values = [value for value in values if value]
    if not observed_values:
        kind = "empty"
    elif _is_numeric(observed_values):
        kind = "numeric"
    else:
        kind = "categorical"
    return TraitColumnSummary(
        name=column,
        kind=kind,
        missing_count=sum(1 for value in values if not value),
        distinct_value_count=len(set(observed_values)),
    )


def validate_traits_table(path: Path, *, taxon_column: str | None = None) -> TraitValidationReport:
    """Validate a trait table and infer deterministic column kinds."""
    table = load_taxon_table(path, taxon_column=taxon_column)
    trait_columns = [
        _summarize_trait_column(table, column)
        for column in table.columns
        if column != table.taxon_column
    ]
    return TraitValidationReport(
        path=table.path,
        format=table.format,
        row_count=table.row_count,
        taxon_column=table.taxon_column,
        trait_columns=trait_columns,
    )
