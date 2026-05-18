from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
from bijux_phylogenetics.reports import (
    SupplementaryAncestralStateRow,
    SupplementaryAncestralStateTableResult,
    write_supplementary_ancestral_state_table,
)


def test_supplementary_ancestral_state_table_surfaces_export_publicly() -> None:
    assert reports_api.SupplementaryAncestralStateRow is SupplementaryAncestralStateRow
    assert (
        reports_api.SupplementaryAncestralStateTableResult
        is SupplementaryAncestralStateTableResult
    )
    assert (
        reports_api.write_supplementary_ancestral_state_table
        is write_supplementary_ancestral_state_table
    )
