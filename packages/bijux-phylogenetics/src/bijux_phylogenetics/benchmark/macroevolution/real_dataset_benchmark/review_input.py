from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.datasets.central_european_seashore_flora import (
    CentralEuropeanSeashoreFloraDataset,
)
from bijux_phylogenetics.datasets.study_inputs import (
    load_taxon_table,
)

from .shared import (
    CONTINUOUS_MISSING_VALUE_TAXON,
    DISCRETE_MISSING_VALUE_TAXON,
    EXTRA_TRAIT_TAXON,
    REMOVED_TREE_TAXON,
)


def write_alignment_review_traits_table(
    path: Path,
    dataset: CentralEuropeanSeashoreFloraDataset,
) -> Path:
    table = load_taxon_table(dataset.traits_path, taxon_column=dataset.taxon_column)
    rows = [
        dict(row)
        for row in table.rows
        if row[dataset.taxon_column] != REMOVED_TREE_TAXON
    ]
    extra_row = dict(rows[0])
    extra_row[dataset.taxon_column] = EXTRA_TRAIT_TAXON
    rows.append(extra_row)
    for row in rows:
        if row[dataset.taxon_column] == CONTINUOUS_MISSING_VALUE_TAXON:
            row[dataset.workflow_continuous_trait] = ""
        if row[dataset.taxon_column] == DISCRETE_MISSING_VALUE_TAXON:
            row[dataset.workflow_discrete_trait] = ""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=table.columns)
        writer.writeheader()
        writer.writerows(rows)
    return path
