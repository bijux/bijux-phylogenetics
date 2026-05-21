from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
from bijux_phylogenetics.reports import (
    SupplementaryTaxonTableResult,
    SupplementaryTaxonTableRow,
    write_supplementary_taxon_table,
)


def test_supplementary_taxon_table_surfaces_export_publicly() -> None:
    assert reports_api.SupplementaryTaxonTableResult is SupplementaryTaxonTableResult
    assert reports_api.SupplementaryTaxonTableRow is SupplementaryTaxonTableRow
    assert (
        reports_api.write_supplementary_taxon_table is write_supplementary_taxon_table
    )
