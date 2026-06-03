from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
from bijux_phylogenetics.reports import (
    SupplementaryCladeSupportRow,
    SupplementaryCladeSupportTableResult,
    write_supplementary_clade_support_table,
)


def test_supplementary_clade_support_table_surfaces_export_publicly() -> None:
    assert reports_api.SupplementaryCladeSupportRow is SupplementaryCladeSupportRow
    assert (
        reports_api.SupplementaryCladeSupportTableResult
        is SupplementaryCladeSupportTableResult
    )
    assert (
        reports_api.write_supplementary_clade_support_table
        is write_supplementary_clade_support_table
    )
