from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

from .contracts import (
    AncestralContinuousChangeBranchRow,
    AncestralContinuousChangeCountRow,
)


def table_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def write_continuous_change_count_table(
    path: Path,
    rows: list[AncestralContinuousChangeCountRow],
) -> Path:
    fieldnames = list(asdict(rows[0]).keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=table_delimiter(path),
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
    return path


def write_continuous_change_branch_table(
    path: Path,
    rows: list[AncestralContinuousChangeBranchRow],
) -> Path:
    fieldnames = list(asdict(rows[0]).keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=table_delimiter(path),
        )
        writer.writeheader()
        for row in rows:
            payload = asdict(row)
            payload["child_descendant_taxa"] = "|".join(row.child_descendant_taxa)
            payload["branch_length"] = (
                "" if row.branch_length is None else row.branch_length
            )
            writer.writerow(payload)
    return path
