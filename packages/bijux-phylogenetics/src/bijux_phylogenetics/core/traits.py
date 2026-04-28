from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class TabularSummary:
    """Summary of a TSV metadata or trait table."""

    path: Path
    row_count: int
    columns: list[str]
    index_column: str
    indexed_values: set[str]


def load_tsv_summary(path: Path) -> TabularSummary:
    """Load a TSV file and expose the first linking column."""
    if not path.exists():
        raise FileNotFoundError(f"table file not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"table has no header row: {path}")
        columns = [column.strip() for column in reader.fieldnames]
        rows = list(reader)

    preferred_columns = ("taxon", "taxa", "tip", "sample", "name", columns[0])
    index_column = next(column for column in preferred_columns if column in columns)
    indexed_values = {str(row.get(index_column, "")).strip() for row in rows if str(row.get(index_column, "")).strip()}
    return TabularSummary(
        path=path,
        row_count=len(rows),
        columns=columns,
        index_column=index_column,
        indexed_values=indexed_values,
    )

