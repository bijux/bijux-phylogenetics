from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias


@dataclass(slots=True)
class TaxonTable:
    """Normalized metadata-style table keyed by a required taxon column."""

    path: Path
    format: str
    row_count: int
    columns: list[str]
    taxon_column: str
    rows: list[dict[str, str]]
    taxa: list[str]

    @property
    def indexed_values(self) -> set[str]:
        return set(self.taxa)

    @property
    def index_column(self) -> str:
        return self.taxon_column


@dataclass(slots=True)
class MetadataColumnCompleteness:
    """Completeness summary for one metadata column."""

    name: str
    missing_count: int
    completeness_fraction: float


@dataclass(slots=True)
class MetadataInspectionReport:
    """Stable summary of a metadata table keyed by taxon."""

    path: Path
    format: str
    row_count: int
    column_count: int
    columns: list[str]
    taxon_column: str
    taxa: list[str]
    column_completeness: list[MetadataColumnCompleteness]


@dataclass(slots=True)
class MetadataJoinRow:
    """One tree tip annotated against a metadata row when available."""

    taxon: str
    matched: bool
    values: dict[str, str]


@dataclass(slots=True)
class MetadataJoinReport:
    """Explicit tip-by-tip join of a metadata table onto a tree."""

    path: Path
    tree_taxa: int
    metadata_rows: int
    taxon_column: str
    joined_rows: list[MetadataJoinRow]
    missing_from_metadata: list[str]
    extra_metadata_taxa: list[str]


@dataclass(slots=True)
class TaxonTableIndexAudit:
    """Index-level audit of a taxon-keyed table before strict loading."""

    path: Path
    format: str
    row_count: int
    taxon_column: str
    taxa: list[str]
    duplicate_taxa: list[str]
    empty_taxon_rows: list[int]


TableValue: TypeAlias = str | int | float | bool | None


__all__ = [
    "MetadataColumnCompleteness",
    "MetadataInspectionReport",
    "MetadataJoinReport",
    "MetadataJoinRow",
    "TableValue",
    "TaxonTable",
    "TaxonTableIndexAudit",
]
