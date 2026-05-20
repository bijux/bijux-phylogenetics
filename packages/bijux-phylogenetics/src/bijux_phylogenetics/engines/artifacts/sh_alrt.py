from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .support import ShAlrtSupportSummaryReport

__all__ = [
    "ShAlrtSupportRow",
    "build_conflicting_sh_alrt_support_rows",
    "build_sh_alrt_support_rows",
    "write_sh_alrt_support_table",
]


@dataclass(frozen=True, slots=True)
class ShAlrtSupportRow:
    """One combined SH-aLRT and UFBoot support record."""

    node: str
    descendant_taxa: tuple[str, ...]
    sh_alrt_support: float | None
    sh_alrt_support_fraction: float | None
    ufboot_support: float | None
    ufboot_support_fraction: float | None
    is_backbone: bool
    sh_alrt_strong: bool
    ufboot_strong: bool
    conflicting_support_signal: bool
    support_agreement: str


def build_sh_alrt_support_rows(
    summary: ShAlrtSupportSummaryReport,
) -> list[ShAlrtSupportRow]:
    """Convert one combined support summary into durable branch rows."""
    return [
        ShAlrtSupportRow(
            node=node.node,
            descendant_taxa=tuple(node.descendant_taxa),
            sh_alrt_support=node.sh_alrt_support,
            sh_alrt_support_fraction=node.sh_alrt_support_fraction,
            ufboot_support=node.ufboot_support,
            ufboot_support_fraction=node.ufboot_support_fraction,
            is_backbone=node.is_backbone,
            sh_alrt_strong=node.sh_alrt_strong,
            ufboot_strong=node.ufboot_strong,
            conflicting_support_signal=node.conflicting_support_signal,
            support_agreement=node.support_agreement,
        )
        for node in summary.nodes
    ]


def build_conflicting_sh_alrt_support_rows(
    summary: ShAlrtSupportSummaryReport,
) -> list[ShAlrtSupportRow]:
    """Return only the branches with conflicting support signals."""
    return [
        row
        for row in build_sh_alrt_support_rows(summary)
        if row.conflicting_support_signal
    ]


def _serialize_taxa(taxa: tuple[str, ...]) -> str:
    return ",".join(taxa)


def _render_numeric(value: float | None) -> str:
    return "" if value is None else format(value, ".12g")


def write_sh_alrt_support_table(
    path: Path,
    rows: list[ShAlrtSupportRow],
) -> Path:
    """Write one combined SH-aLRT/UFBoot branch-support table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "\t".join(
            [
                "node",
                "descendant_taxa",
                "sh_alrt_support",
                "sh_alrt_support_fraction",
                "ufboot_support",
                "ufboot_support_fraction",
                "is_backbone",
                "sh_alrt_strong",
                "ufboot_strong",
                "conflicting_support_signal",
                "support_agreement",
            ]
        )
    ]
    lines.extend(
        "\t".join(
            [
                row.node,
                _serialize_taxa(row.descendant_taxa),
                _render_numeric(row.sh_alrt_support),
                _render_numeric(row.sh_alrt_support_fraction),
                _render_numeric(row.ufboot_support),
                _render_numeric(row.ufboot_support_fraction),
                "true" if row.is_backbone else "false",
                "true" if row.sh_alrt_strong else "false",
                "true" if row.ufboot_strong else "false",
                "true" if row.conflicting_support_signal else "false",
                row.support_agreement,
            ]
        )
        for row in rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
