from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .support import BootstrapSupportSummaryReport

__all__ = [
    "BootstrapSupportHistogramRow",
    "BootstrapSupportRow",
    "build_bootstrap_support_histogram_rows",
    "build_bootstrap_support_rows",
    "build_low_support_bootstrap_rows",
    "classify_bootstrap_support_bucket",
    "write_bootstrap_support_histogram",
    "write_bootstrap_support_table",
]


@dataclass(frozen=True, slots=True)
class BootstrapSupportRow:
    """One branch-level support record for one bootstrap-supported tree."""

    node: str
    descendant_taxa: tuple[str, ...]
    support: float
    support_fraction: float
    is_backbone: bool
    support_bucket: str
    low_support: bool


@dataclass(frozen=True, slots=True)
class BootstrapSupportHistogramRow:
    """One support-frequency bucket for reviewer-facing histogram export."""

    support_bucket: str
    minimum_support: float | None
    maximum_support: float | None
    node_count: int


def classify_bootstrap_support_bucket(support: float) -> str:
    """Map one numeric support value to the governed histogram bucket."""
    if support < 50.0:
        return "lt50"
    if support < 70.0:
        return "50to69"
    if support < 90.0:
        return "70to89"
    return "ge90"


def build_bootstrap_support_rows(
    summary: BootstrapSupportSummaryReport,
    *,
    low_support_threshold: float = 70.0,
) -> list[BootstrapSupportRow]:
    """Convert one bootstrap support summary into durable branch rows."""
    return [
        BootstrapSupportRow(
            node=node.node,
            descendant_taxa=tuple(node.descendant_taxa),
            support=node.support,
            support_fraction=node.support_fraction,
            is_backbone=node.is_backbone,
            support_bucket=classify_bootstrap_support_bucket(node.support),
            low_support=node.support < low_support_threshold,
        )
        for node in summary.nodes
    ]


def build_low_support_bootstrap_rows(
    summary: BootstrapSupportSummaryReport,
    *,
    low_support_threshold: float = 70.0,
) -> list[BootstrapSupportRow]:
    """Return only the weakly supported clades from one support summary."""
    return [
        row
        for row in build_bootstrap_support_rows(
            summary,
            low_support_threshold=low_support_threshold,
        )
        if row.low_support
    ]


def build_bootstrap_support_histogram_rows(
    summary: BootstrapSupportSummaryReport,
) -> list[BootstrapSupportHistogramRow]:
    """Convert one support histogram into stable reviewer-facing rows."""
    return [
        BootstrapSupportHistogramRow(
            support_bucket="lt50",
            minimum_support=None,
            maximum_support=49.999999999999,
            node_count=summary.support_histogram["lt50"],
        ),
        BootstrapSupportHistogramRow(
            support_bucket="50to69",
            minimum_support=50.0,
            maximum_support=69.999999999999,
            node_count=summary.support_histogram["50to69"],
        ),
        BootstrapSupportHistogramRow(
            support_bucket="70to89",
            minimum_support=70.0,
            maximum_support=89.999999999999,
            node_count=summary.support_histogram["70to89"],
        ),
        BootstrapSupportHistogramRow(
            support_bucket="ge90",
            minimum_support=90.0,
            maximum_support=None,
            node_count=summary.support_histogram["ge90"],
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


def write_bootstrap_support_table(
    path: Path,
    rows: list[BootstrapSupportRow],
) -> Path:
    """Write branch-level support evidence for one bootstrap workflow."""
    return _write_tsv(
        path,
        header=[
            "node",
            "descendant_taxa",
            "support",
            "support_fraction",
            "is_backbone",
            "support_bucket",
            "low_support",
        ],
        rows=[
            [
                row.node,
                _serialize_support_taxa(row.descendant_taxa),
                format(row.support, ".12g"),
                format(row.support_fraction, ".12g"),
                "true" if row.is_backbone else "false",
                row.support_bucket,
                "true" if row.low_support else "false",
            ]
            for row in rows
        ],
    )


def write_bootstrap_support_histogram(
    path: Path,
    rows: list[BootstrapSupportHistogramRow],
) -> Path:
    """Write the support histogram used for reviewer-facing bootstrap review."""
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
