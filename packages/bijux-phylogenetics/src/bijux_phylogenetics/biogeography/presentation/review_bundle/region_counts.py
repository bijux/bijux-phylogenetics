from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree

from .contracts import BiogeographyRegionCountRow


def summarize_biogeography_region_counts(
    tree_path: Path,
    table_path: Path,
    *,
    trait: str,
    taxon_column: str,
    excluded_taxa: set[str],
) -> list[BiogeographyRegionCountRow]:
    """Count observed regions after tree overlap and exclusion auditing."""
    tree = load_tree(tree_path)
    table = load_taxon_table(table_path, taxon_column=taxon_column)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    counts: dict[str, int] = {}
    analyzed_taxa = 0
    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None or taxon in excluded_taxa:
            continue
        region = row.get(trait, "").strip()
        if not region:
            continue
        counts[region] = counts.get(region, 0) + 1
        analyzed_taxa += 1
    if analyzed_taxa == 0:
        return []
    return [
        BiogeographyRegionCountRow(
            region=region,
            tip_taxon_count=count,
            analyzed_taxon_fraction=count / analyzed_taxa,
        )
        for region, count in sorted(counts.items())
    ]


def write_biogeography_region_count_table(
    path: Path,
    rows: list[BiogeographyRegionCountRow],
) -> Path:
    """Write one observed-region count ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "region",
            "tip_taxon_count",
            "analyzed_taxon_fraction",
        ],
        rows=[
            {
                "region": row.region,
                "tip_taxon_count": str(row.tip_taxon_count),
                "analyzed_taxon_fraction": str(row.analyzed_taxon_fraction),
            }
            for row in rows
        ],
    )
