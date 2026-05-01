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


@dataclass(slots=True)
class TraitTablePruningReport:
    """Explicit record of pruning a trait table to tree taxa."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    original_row_count: int
    kept_taxa: list[str]
    removed_taxa: list[str]


@dataclass(slots=True)
class MissingTraitValue:
    """One missing trait observation tied to a taxon and column."""

    taxon: str
    trait: str


@dataclass(slots=True)
class TraitMissingValueReport:
    """Explicit missing trait calls by taxon and trait column."""

    path: Path
    taxon_column: str
    missing_values: list[MissingTraitValue]


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


def _is_binary(values: list[str]) -> bool:
    normalized = {value.strip().lower() for value in values}
    binary_tokens = {
        "0",
        "1",
        "false",
        "true",
        "no",
        "yes",
        "absent",
        "present",
    }
    return bool(normalized) and normalized <= binary_tokens and len(normalized) <= 2


def _is_text(values: list[str], *, row_count: int) -> bool:
    distinct_values = len(set(values))
    if any(" " in value for value in values):
        return True
    return distinct_values > max(3, row_count // 2)


def _summarize_trait_column(table: TaxonTable, column: str) -> TraitColumnSummary:
    values = [row[column] for row in table.rows]
    observed_values = [value for value in values if value]
    if not observed_values:
        kind = "empty"
    elif _is_numeric(observed_values):
        kind = "numeric"
    elif _is_binary(observed_values):
        kind = "binary"
    elif _is_text(observed_values, row_count=len(table.rows)):
        kind = "text"
    else:
        kind = "categorical"
    return TraitColumnSummary(
        name=column,
        kind=kind,
        missing_count=sum(1 for value in values if not value),
        missing_fraction=sum(1 for value in values if not value) / max(len(values), 1),
        distinct_value_count=len(set(observed_values)),
    )


def validate_traits_table(
    path: Path, *, taxon_column: str | None = None
) -> TraitValidationReport:
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


def prune_traits_to_tree(
    tree_path: Path,
    traits_path: Path,
    *,
    taxon_column: str | None = None,
) -> tuple[list[dict[str, str]], TraitTablePruningReport]:
    """Prune a trait table to the taxa present in a tree while preserving tree tip order."""
    tree = load_tree(tree_path)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    kept_rows = [
        dict(rows_by_taxon[taxon]) for taxon in tree.tip_names if taxon in rows_by_taxon
    ]
    if not kept_rows:
        raise MetadataJoinError(
            "no overlapping taxa remain after trait pruning request"
        )

    kept_taxa = [row[table.taxon_column] for row in kept_rows]
    removed_taxa = sorted(set(table.taxa) - set(kept_taxa))
    return kept_rows, TraitTablePruningReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        original_row_count=table.row_count,
        kept_taxa=kept_taxa,
        removed_taxa=removed_taxa,
    )


def detect_missing_trait_values(
    path: Path,
    *,
    taxon_column: str | None = None,
) -> TraitMissingValueReport:
    """Return every missing trait value with its taxon and column name."""
    table = load_taxon_table(path, taxon_column=taxon_column)
    missing_values = [
        MissingTraitValue(
            taxon=row[table.taxon_column],
            trait=column,
        )
        for row in table.rows
        for column in table.columns
        if column != table.taxon_column and not row[column]
    ]
    return TraitMissingValueReport(
        path=table.path,
        taxon_column=table.taxon_column,
        missing_values=missing_values,
    )
