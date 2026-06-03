from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
from bijux_phylogenetics.reports import (
    SupplementaryDiversificationRow,
    SupplementaryDiversificationTableResult,
    write_supplementary_diversification_table,
)


def test_supplementary_diversification_table_surfaces_export_publicly() -> None:
    assert (
        reports_api.SupplementaryDiversificationRow is SupplementaryDiversificationRow
    )
    assert (
        reports_api.SupplementaryDiversificationTableResult
        is SupplementaryDiversificationTableResult
    )
    assert (
        reports_api.write_supplementary_diversification_table
        is write_supplementary_diversification_table
    )
