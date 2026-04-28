from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.metadata import TaxonTable, load_taxon_table
from bijux_phylogenetics.errors import MetadataJoinError
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class TraitColumnSummary:
    """Deterministic schema summary for one trait column."""

    name: str
    kind: str
    missing_count: int
    missing_fraction: float
    distinct_value_count: int


@dataclass(slots=True)
class TraitValidationReport:
    """Stable summary of a validated trait table."""

    path: Path
    format: str
    row_count: int
    taxon_column: str
    trait_columns: list[TraitColumnSummary]


@dataclass(slots=True)
class TraitLinkageReport:
    """Summary of how a trait table joins against a tree tip set."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    tree_taxa: int
    trait_taxa: int
    linked_taxa: int
    usable_taxa: list[str]
    missing_from_traits: list[str]
    extra_trait_taxa: list[str]


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
        missing_fraction=sum(1 for value in values if not value) / max(len(values), 1),
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


def link_tree_to_traits(
    tree_path: Path,
    traits_path: Path,
    *,
    taxon_column: str | None = None,
    strict: bool = False,
) -> TraitLinkageReport:
    """Report how a traits table links against tree tips."""
    tree = load_tree(tree_path)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    tree_taxa = set(tree.tip_names)
    trait_taxa = set(table.taxa)
    missing_from_traits = sorted(tree_taxa - trait_taxa)
    extra_trait_taxa = sorted(trait_taxa - tree_taxa)

    if strict and (missing_from_traits or extra_trait_taxa):
        raise MetadataJoinError(
            "trait linkage mismatch: "
            f"{len(missing_from_traits)} tree taxa missing from traits and "
            f"{len(extra_trait_taxa)} trait taxa absent from tree"
        )

    usable_taxa = sorted(tree_taxa & trait_taxa)
    return TraitLinkageReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        tree_taxa=len(tree_taxa),
        trait_taxa=len(trait_taxa),
        linked_taxa=len(usable_taxa),
        usable_taxa=usable_taxa,
        missing_from_traits=missing_from_traits,
        extra_trait_taxa=extra_trait_taxa,
    )


def detect_unusable_trait_columns(
    path: Path,
    *,
    missingness_threshold: float,
    taxon_column: str | None = None,
) -> list[TraitColumnSummary]:
    """Return trait columns whose missingness exceeds the given threshold."""
    if not 0.0 <= missingness_threshold <= 1.0:
        raise ValueError(
            f"missingness threshold must be between 0 and 1 inclusive, got {missingness_threshold}"
        )
    report = validate_traits_table(path, taxon_column=taxon_column)
    return [
        column
        for column in report.trait_columns
        if column.missing_fraction > missingness_threshold
    ]
