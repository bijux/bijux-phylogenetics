from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import TaxonTable


def row_lookup(table: TaxonTable) -> dict[str, dict[str, str]]:
    return {row[table.taxon_column]: row for row in table.rows}


def stringify_list(values: list[str]) -> str:
    return "|".join(values)


def stringify_mapping(values: dict[str, str]) -> str:
    return "|".join(f"{key}={value}" for key, value in sorted(values.items()))


def table_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def build_dynamic_columns(
    metadata_table: TaxonTable,
    traits_table: TaxonTable,
) -> tuple[list[str], list[str]]:
    metadata_columns = [
        f"metadata_{column}"
        for column in metadata_table.columns
        if column != metadata_table.taxon_column
    ]
    trait_columns = [
        f"trait_{column}"
        for column in traits_table.columns
        if column != traits_table.taxon_column
    ]
    return metadata_columns, trait_columns


def write_dict_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[dict[str, object | str]],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter=table_delimiter(path),
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path
