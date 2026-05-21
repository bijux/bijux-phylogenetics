from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

from .contracts import TreeBranchStatisticsRow, TreeSupportRow


def table_delimiter(path: Path) -> str:
    """Return the governed delimiter for one artifact table path."""
    return "," if path.suffix.lower() == ".csv" else "\t"


def write_tree_support_table(path: Path, rows: list[TreeSupportRow]) -> Path:
    """Write the reviewer-facing support table artifact."""
    fieldnames = [
        "node_kind",
        "node",
        "node_label",
        "descendant_taxa",
        "support",
        "support_fraction",
        "support_class",
        "branch_length",
        "root_depth",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter=table_delimiter(path)
        )
        writer.writeheader()
        for row in rows:
            payload = asdict(row)
            payload["node_label"] = "" if row.node_label is None else row.node_label
            payload["descendant_taxa"] = "|".join(row.descendant_taxa)
            payload["support"] = "" if row.support is None else row.support
            payload["support_fraction"] = (
                "" if row.support_fraction is None else row.support_fraction
            )
            payload["branch_length"] = (
                "" if row.branch_length is None else row.branch_length
            )
            payload["root_depth"] = "" if row.root_depth is None else row.root_depth
            writer.writerow(payload)
    return path


def write_tree_branch_statistics_table(
    path: Path, row: TreeBranchStatisticsRow
) -> Path:
    """Write the reviewer-facing branch statistics artifact."""
    fieldnames = list(asdict(row).keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter=table_delimiter(path)
        )
        writer.writeheader()
        writer.writerow(asdict(row))
    return path
