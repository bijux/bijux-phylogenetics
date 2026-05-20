from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .support import FastTreeSupportSummaryReport

__all__ = [
    "FastTreeSupportHistogramRow",
    "FastTreeSupportRow",
    "build_fasttree_low_support_rows",
    "build_fasttree_support_histogram_rows",
    "build_fasttree_support_rows",
    "classify_fasttree_support_bucket",
    "write_fasttree_support_histogram",
    "write_fasttree_support_table",
]


@dataclass(frozen=True, slots=True)
class FastTreeSupportRow:
    """One branch-level FastTree local-support record."""

    node: str
    descendant_taxa: tuple[str, ...]
    local_support: float
    support_fraction: float
    is_backbone: bool
    support_bucket: str
    low_support: bool


@dataclass(frozen=True, slots=True)
class FastTreeSupportHistogramRow:
    """One reviewer-facing bucket for FastTree local support."""

    support_bucket: str
    minimum_support: float | None
    maximum_support: float | None
    node_count: int


def classify_fasttree_support_bucket(local_support: float) -> str:
    """Map one FastTree local-support proportion to a stable review bucket."""
    if local_support < 0.5:
        return "lt0p5"
    if local_support < 0.7:
        return "0p5to0p69"
    if local_support < 0.9:
        return "0p7to0p89"
    return "ge0p9"


def build_fasttree_support_rows(
    summary: FastTreeSupportSummaryReport,
    *,
    low_support_threshold: float = 0.7,
) -> list[FastTreeSupportRow]:
    """Convert one FastTree support summary into durable branch rows."""
    return [
        FastTreeSupportRow(
            node=node.node,
            descendant_taxa=tuple(node.descendant_taxa),
            local_support=node.local_support,
            support_fraction=node.support_fraction,
            is_backbone=node.is_backbone,
            support_bucket=classify_fasttree_support_bucket(node.local_support),
            low_support=node.local_support < low_support_threshold,
        )
        for node in summary.nodes
    ]


def build_fasttree_low_support_rows(
    summary: FastTreeSupportSummaryReport,
    *,
    low_support_threshold: float = 0.7,
) -> list[FastTreeSupportRow]:
    """Return only weakly supported FastTree clades."""
    return [
        row
        for row in build_fasttree_support_rows(
            summary,
            low_support_threshold=low_support_threshold,
        )
        if row.low_support
    ]


def build_fasttree_support_histogram_rows(
    summary: FastTreeSupportSummaryReport,
) -> list[FastTreeSupportHistogramRow]:
    """Convert one FastTree support histogram into stable reviewer-facing rows."""
    return [
        FastTreeSupportHistogramRow(
            support_bucket="lt0p5",
            minimum_support=None,
            maximum_support=0.499999999999,
            node_count=summary.support_histogram["lt0p5"],
        ),
        FastTreeSupportHistogramRow(
            support_bucket="0p5to0p69",
            minimum_support=0.5,
            maximum_support=0.699999999999,
            node_count=summary.support_histogram["0p5to0p69"],
        ),
        FastTreeSupportHistogramRow(
            support_bucket="0p7to0p89",
            minimum_support=0.7,
            maximum_support=0.899999999999,
            node_count=summary.support_histogram["0p7to0p89"],
        ),
        FastTreeSupportHistogramRow(
            support_bucket="ge0p9",
            minimum_support=0.9,
            maximum_support=None,
            node_count=summary.support_histogram["ge0p9"],
        ),
    ]


def _serialize_support_taxa(taxa: tuple[str, ...]) -> str:
    return ",".join(taxa)


def _write_tsv(path: Path, *, header: list[str], rows: list[list[str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_fasttree_support_table(path: Path, rows: list[FastTreeSupportRow]) -> Path:
    """Write branch-level FastTree support evidence for one workflow."""
    return _write_tsv(
        path,
        header=[
            "node",
            "descendant_taxa",
            "local_support",
            "support_fraction",
            "is_backbone",
            "support_bucket",
            "low_support",
        ],
        rows=[
            [
                row.node,
                _serialize_support_taxa(row.descendant_taxa),
                format(row.local_support, ".12g"),
                format(row.support_fraction, ".12g"),
                "true" if row.is_backbone else "false",
                row.support_bucket,
                "true" if row.low_support else "false",
            ]
            for row in rows
        ],
    )


def write_fasttree_support_histogram(
    path: Path,
    rows: list[FastTreeSupportHistogramRow],
) -> Path:
    """Write the reviewer-facing FastTree local-support histogram."""
    return _write_tsv(
        path,
        header=[
            "support_bucket",
            "minimum_support",
            "maximum_support",
            "node_count",
        ],
        rows=[
            [
                row.support_bucket,
                ""
                if row.minimum_support is None
                else format(row.minimum_support, ".12g"),
                ""
                if row.maximum_support is None
                else format(row.maximum_support, ".12g"),
                str(row.node_count),
            ]
            for row in rows
        ],
    )
