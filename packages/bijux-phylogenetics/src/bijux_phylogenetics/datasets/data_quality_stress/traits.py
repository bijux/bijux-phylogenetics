from __future__ import annotations

import csv
from pathlib import Path

from .models import TraitDuplicateResolution, TraitMissingObservation


def load_permissive_trait_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"trait table has no header row: {path}")
        columns = [column.strip() for column in reader.fieldnames]
        return [
            {column: str(row.get(column, "")).strip() for column in columns}
            for row in reader
        ]


def selected_trait_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[tuple[int, dict[str, str]]]] = {}
    for index, row in enumerate(rows, start=2):
        grouped.setdefault(row["taxon"], []).append((index, row))
    selected: list[tuple[int, dict[str, str]]] = []
    for entries in grouped.values():
        selected.append(select_trait_row(entries))
    return [dict(row) for _, row in sorted(selected, key=lambda item: item[0])]


def resolve_duplicate_traits(
    rows: list[dict[str, str]],
) -> list[TraitDuplicateResolution]:
    grouped: dict[str, list[tuple[int, dict[str, str]]]] = {}
    for index, row in enumerate(rows, start=2):
        grouped.setdefault(row["taxon"], []).append((index, row))
    resolutions: list[TraitDuplicateResolution] = []
    for taxon, entries in sorted(grouped.items()):
        if len(entries) < 2:
            continue
        selected_row_number, selected_row = select_trait_row(entries)
        resolutions.append(
            TraitDuplicateResolution(
                taxon=taxon,
                occurrence_count=len(entries),
                selected_row_number=selected_row_number,
                selected_non_missing_field_count=non_missing_field_count(selected_row),
                discarded_row_numbers=[
                    row_number
                    for row_number, _ in entries
                    if row_number != selected_row_number
                ],
                selected_reason="highest_non_missing_field_count_then_first_row",
            )
        )
    return resolutions


def select_trait_row(
    entries: list[tuple[int, dict[str, str]]],
) -> tuple[int, dict[str, str]]:
    return max(
        entries,
        key=lambda item: (non_missing_field_count(item[1]), -item[0]),
    )


def non_missing_field_count(row: dict[str, str]) -> int:
    return sum(1 for key, value in row.items() if key != "taxon" and value)


def detect_missing_traits(
    rows: list[dict[str, str]],
    *,
    required_traits: set[str],
    duplicate_lookup: dict[str, TraitDuplicateResolution],
) -> list[TraitMissingObservation]:
    observations: list[TraitMissingObservation] = []
    for row_number, row in enumerate(rows, start=2):
        for trait, value in row.items():
            if trait == "taxon" or value:
                continue
            required_for_analysis = trait in required_traits
            if row["taxon"] in duplicate_lookup and (
                row_number != duplicate_lookup[row["taxon"]].selected_row_number
            ):
                action = "dropped_duplicate_row"
            elif required_for_analysis:
                action = "drop_taxon_from_cleaned_traits"
            else:
                action = "preserve_nonrequired_missingness"
            observations.append(
                TraitMissingObservation(
                    taxon=row["taxon"],
                    row_number=row_number,
                    trait=trait,
                    required_for_analysis=required_for_analysis,
                    action=action,
                )
            )
    return observations
