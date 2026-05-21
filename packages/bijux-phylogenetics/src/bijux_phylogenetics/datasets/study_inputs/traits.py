from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import MetadataJoinError

from .models import TableValue, TaxonTable
from .tables import load_taxon_table, write_taxon_rows


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


@dataclass(slots=True)
class TaxonNameMismatchRow:
    """One taxon-name mismatch between a tree tip set and a trait table."""

    mismatch_side: str
    taxon: str
    present_in_tree: bool
    present_in_traits: bool


@dataclass(slots=True)
class TreeTraitNameCheckReport:
    """Stable taxon-name mismatch report aligned to `geiger::name.check`."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    tree_taxa: int
    trait_taxa: int
    tree_not_data: list[str]
    data_not_tree: list[str]
    mismatch_rows: list[TaxonNameMismatchRow]
    compatible: bool
    reference_outcome: str
    matching_policy: str


@dataclass(slots=True)
class TreeTraitAlignmentReport:
    """Stable record of aligning one tree to one taxon-keyed trait table."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    original_tree_taxa: int
    original_trait_taxa: int
    aligned_taxa: list[str]
    dropped_tree_taxa: list[str]
    dropped_trait_taxa: list[str]
    dropped_missing_value_taxa: list[str]
    missing_value_calls: list[MissingTraitValue]
    tree_drop_policy: str
    trait_drop_policy: str
    missing_value_policy: str


@dataclass(slots=True)
class TreeTraitAlignment:
    """Pruned tree plus tree-ordered trait rows after explicit taxon alignment."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    tree: PhyloTree
    rows: list[dict[str, str]]
    report: TreeTraitAlignmentReport


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


def check_tree_and_trait_taxon_names(
    tree_path: Path,
    traits_path: Path,
    *,
    taxon_column: str | None = None,
) -> TreeTraitNameCheckReport:
    """Report tree-versus-trait taxon mismatches with `geiger::name.check` semantics."""
    tree = load_tree(tree_path)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    trait_taxa = set(table.taxa)
    tree_not_data = [taxon for taxon in tree.tip_names if taxon not in trait_taxa]
    data_not_tree = [taxon for taxon in table.taxa if taxon not in set(tree.tip_names)]
    mismatch_rows = [
        TaxonNameMismatchRow(
            mismatch_side="tree_not_data",
            taxon=taxon,
            present_in_tree=True,
            present_in_traits=False,
        )
        for taxon in tree_not_data
    ] + [
        TaxonNameMismatchRow(
            mismatch_side="data_not_tree",
            taxon=taxon,
            present_in_tree=False,
            present_in_traits=True,
        )
        for taxon in data_not_tree
    ]
    compatible = not mismatch_rows
    return TreeTraitNameCheckReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        tree_taxa=len(tree.tip_names),
        trait_taxa=table.row_count,
        tree_not_data=tree_not_data,
        data_not_tree=data_not_tree,
        mismatch_rows=mismatch_rows,
        compatible=compatible,
        reference_outcome="OK" if compatible else "mismatch",
        matching_policy="case-sensitive-exact-label-matching",
    )


def write_tree_trait_name_mismatch_table(
    path: Path,
    report: TreeTraitNameCheckReport,
) -> Path:
    """Write one machine-readable taxon-name mismatch table."""
    rows: list[dict[str, TableValue]] = [
        {
            "mismatch_side": row.mismatch_side,
            "taxon": row.taxon,
            "present_in_tree": row.present_in_tree,
            "present_in_traits": row.present_in_traits,
        }
        for row in report.mismatch_rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "mismatch_side",
            "taxon",
            "present_in_tree",
            "present_in_traits",
        ],
        rows=rows,
    )


def link_tree_to_traits(
    tree_path: Path,
    traits_path: Path,
    *,
    taxon_column: str | None = None,
    strict: bool = False,
) -> TraitLinkageReport:
    """Report how a traits table links against tree tips."""
    alignment = align_tree_and_trait_table(
        tree_path,
        traits_path,
        taxon_column=taxon_column,
    )
    report = alignment.report
    missing_from_traits = report.dropped_tree_taxa
    extra_trait_taxa = report.dropped_trait_taxa

    if strict and (missing_from_traits or extra_trait_taxa):
        tree_taxa = sorted(set(report.aligned_taxa) | set(report.dropped_tree_taxa))
        trait_taxa = sorted(set(report.aligned_taxa) | set(report.dropped_trait_taxa))
        raise MetadataJoinError(
            "trait linkage mismatch: "
            f"{len(missing_from_traits)} tree taxa missing from traits and "
            f"{len(extra_trait_taxa)} trait taxa absent from tree",
            details={
                "tree_path": str(tree_path),
                "traits_path": str(traits_path),
                "taxon_column": report.taxon_column,
                "missing_from_traits": missing_from_traits,
                "extra_trait_taxa": extra_trait_taxa,
                "failure_reason": "tree_trait_taxon_mismatch",
                "scientific_explanation": (
                    "The tree and trait table do not describe the same biological taxon set, "
                    "so comparative interpretation would silently drop or misalign taxa."
                ),
                "likely_causes": [
                    "one or more tree tips are missing from the trait table",
                    "the trait table contains taxa that are absent from the tree",
                ],
                "actionable_fixes": [
                    "add the missing tree taxa to the trait table or prune the tree intentionally",
                    "remove extra trait-table taxa that are not represented in the tree",
                ],
                "evidence": {
                    "tree_taxa": sorted(tree_taxa),
                    "trait_taxa": sorted(trait_taxa),
                    "missing_from_traits": missing_from_traits,
                    "extra_trait_taxa": extra_trait_taxa,
                },
            },
        )

    return TraitLinkageReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=report.taxon_column,
        tree_taxa=report.original_tree_taxa,
        trait_taxa=report.original_trait_taxa,
        linked_taxa=len(report.aligned_taxa),
        usable_taxa=list(report.aligned_taxa),
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
    alignment = align_tree_and_trait_table(
        tree_path,
        traits_path,
        taxon_column=taxon_column,
    )
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    kept_rows = [dict(row) for row in alignment.rows]
    kept_taxa = list(alignment.report.aligned_taxa)
    removed_taxa = list(alignment.report.dropped_trait_taxa)
    return kept_rows, TraitTablePruningReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        original_row_count=table.row_count,
        kept_taxa=kept_taxa,
        removed_taxa=removed_taxa,
    )


def align_tree_and_trait_table(
    tree_path: Path,
    traits_path: Path,
    *,
    taxon_column: str | None = None,
    required_trait_columns: Sequence[str] = (),
    drop_missing_for_columns: Sequence[str] = (),
) -> TreeTraitAlignment:
    """Align a tree and trait table over shared taxa while preserving tree tip order."""
    tree = load_tree(tree_path)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    for column in required_trait_columns:
        if column not in table.columns:
            raise MetadataJoinError(f"trait table does not contain column '{column}'")
    for column in drop_missing_for_columns:
        if column not in table.columns:
            raise MetadataJoinError(f"trait table does not contain column '{column}'")

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    tree_taxa = set(tree.tip_names)
    trait_taxa = set(table.taxa)
    dropped_tree_taxa = sorted(tree_taxa - trait_taxa)
    dropped_trait_taxa = sorted(trait_taxa - tree_taxa)

    missing_value_columns = tuple(
        column for column in drop_missing_for_columns if column != table.taxon_column
    )
    missing_value_calls: list[MissingTraitValue] = []
    dropped_missing_value_taxa: list[str] = []
    aligned_rows: list[dict[str, str]] = []
    aligned_taxa: list[str] = []

    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None:
            continue
        missing_columns = [
            column for column in missing_value_columns if not row[column]
        ]
        if missing_columns:
            dropped_missing_value_taxa.append(taxon)
            missing_value_calls.extend(
                MissingTraitValue(taxon=taxon, trait=column)
                for column in missing_columns
            )
            continue
        aligned_taxa.append(taxon)
        aligned_rows.append(dict(row))

    if not aligned_taxa:
        raise MetadataJoinError("no overlapping taxa remain after trait alignment")

    pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, aligned_taxa)
    report = TreeTraitAlignmentReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        original_tree_taxa=len(tree.tip_names),
        original_trait_taxa=table.row_count,
        aligned_taxa=list(pruned_tree.tip_names),
        dropped_tree_taxa=dropped_tree_taxa,
        dropped_trait_taxa=dropped_trait_taxa,
        dropped_missing_value_taxa=sorted(dropped_missing_value_taxa),
        missing_value_calls=missing_value_calls,
        tree_drop_policy="drop-tree-tips-absent-from-traits",
        trait_drop_policy="drop-trait-rows-absent-from-tree",
        missing_value_policy=(
            "retain-overlapping-missing-values"
            if not missing_value_columns
            else "drop-overlapping-missing-values-for-requested-traits"
        ),
    )
    return TreeTraitAlignment(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        tree=pruned_tree,
        rows=aligned_rows,
        report=report,
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
