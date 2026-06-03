from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import MetadataJoinError


def load_tip_dates_for_tree(
    metadata_path: Path,
    *,
    tree_taxa: list[str],
    taxon_column: str | None,
    date_column: str,
) -> tuple[dict[str, float], str]:
    """Load one numeric tip-date table and verify exact taxon coverage."""
    table = load_taxon_table(metadata_path, taxon_column=taxon_column)
    if date_column not in table.columns:
        raise MetadataJoinError(
            f"missing date column '{date_column}' in {metadata_path}"
        )
    tree_taxa_set = set(tree_taxa)
    table_taxa_set = set(table.taxa)
    missing_tree_taxa = sorted(tree_taxa_set - table_taxa_set)
    if missing_tree_taxa:
        raise MetadataJoinError(
            "tip-date table is missing tree taxa: " + ", ".join(missing_tree_taxa)
        )
    extra_table_taxa = sorted(table_taxa_set - tree_taxa_set)
    if extra_table_taxa:
        raise MetadataJoinError(
            "tip-date table contains taxa absent from the tree: "
            + ", ".join(extra_table_taxa)
        )
    tip_dates: dict[str, float] = {}
    for row in table.rows:
        taxon = row[table.taxon_column]
        raw_value = row.get(date_column, "")
        if not raw_value:
            raise MetadataJoinError(
                f"taxon '{taxon}' is missing a numeric date in column '{date_column}'"
            )
        try:
            tip_dates[taxon] = float(raw_value)
        except ValueError as error:
            raise MetadataJoinError(
                f"taxon '{taxon}' has a non-numeric date '{raw_value}'"
            ) from error
    return tip_dates, table.taxon_column


def validate_tip_dates_against_tree(
    tree: PhyloTree,
    tip_dates: Mapping[str, float],
) -> None:
    """Validate that one tip-date mapping covers exactly the tree tip set."""
    tree_taxa = set(tree.tip_names)
    tip_date_taxa = set(tip_dates)
    missing_tree_taxa = sorted(tree_taxa - tip_date_taxa)
    if missing_tree_taxa:
        raise MetadataJoinError(
            "tip-date mapping is missing tree taxa: " + ", ".join(missing_tree_taxa)
        )
    extra_tip_taxa = sorted(tip_date_taxa - tree_taxa)
    if extra_tip_taxa:
        raise MetadataJoinError(
            "tip-date mapping contains taxa absent from the tree: "
            + ", ".join(extra_tip_taxa)
        )
