from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.metadata import MetadataColumnCompleteness, inspect_metadata_table, load_taxon_table
from bijux_phylogenetics.core.traits import detect_unusable_trait_columns, link_tree_to_traits
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class DatasetReadinessSummary:
    """Readiness summary for a tree plus linked metadata and trait tables."""

    tree_path: Path
    metadata_path: Path
    traits_path: Path
    tree_taxa: int
    analysis_taxa: list[str]
    missing_metadata_taxa: list[str]
    missing_trait_taxa: list[str]
    metadata_only_taxa: list[str]
    trait_only_taxa: list[str]
    metadata_column_completeness: list[MetadataColumnCompleteness]
    unusable_trait_columns: list[str]
    ready_for_comparative_analysis: bool
    blockers: list[str]
    warnings: list[str]


def summarize_dataset_readiness(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    trait_missingness_threshold: float = 0.25,
) -> DatasetReadinessSummary:
    """Summarize whether a tree plus linked tables are ready for comparative analysis."""
    tree_validation = validate_tree_path(tree_path, strict=True)
    metadata = inspect_metadata_table(metadata_path)
    tree_taxa = set(load_tree(tree_path).tip_names)
    metadata_table = load_taxon_table(metadata_path)
    metadata_taxa = set(metadata_table.taxa)
    traits_linkage = link_tree_to_traits(tree_path, traits_path)
    unusable_trait_columns = detect_unusable_trait_columns(
        traits_path,
        missingness_threshold=trait_missingness_threshold,
    )

    analysis_taxa = sorted(
        tree_taxa.intersection(metadata_taxa, traits_linkage.usable_taxa)
    )
    blockers: list[str] = []
    warnings: list[str] = []

    if tree_validation.branch_length_status != "complete":
        blockers.append("tree requires complete branch lengths")
    missing_metadata_taxa = sorted(tree_taxa - metadata_taxa)
    metadata_only_taxa = sorted(metadata_taxa - tree_taxa)
    if missing_metadata_taxa:
        blockers.append("metadata table is missing one or more tree taxa")
    if traits_linkage.missing_from_traits:
        blockers.append("trait table is missing one or more tree taxa")
    if unusable_trait_columns:
        blockers.append("one or more trait columns exceed the missingness threshold")
    if len(analysis_taxa) < 2:
        blockers.append("fewer than two taxa remain after intersecting tree, metadata, and traits")
    if metadata_only_taxa:
        warnings.append("metadata table contains taxa absent from the tree")
    if traits_linkage.extra_trait_taxa:
        warnings.append("trait table contains taxa absent from the tree")

    return DatasetReadinessSummary(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        tree_taxa=tree_validation.tip_count,
        analysis_taxa=analysis_taxa,
        missing_metadata_taxa=missing_metadata_taxa,
        missing_trait_taxa=traits_linkage.missing_from_traits,
        metadata_only_taxa=metadata_only_taxa,
        trait_only_taxa=traits_linkage.extra_trait_taxa,
        metadata_column_completeness=metadata.column_completeness,
        unusable_trait_columns=[column.name for column in unusable_trait_columns],
        ready_for_comparative_analysis=not blockers,
        blockers=blockers,
        warnings=warnings,
    )
