from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
from bijux_phylogenetics.reports import (
    SupplementaryComparativeModelRow,
    SupplementaryComparativeModelTableResult,
    write_supplementary_comparative_model_table,
)


def test_supplementary_comparative_model_table_surfaces_export_publicly() -> None:
    assert (
        reports_api.SupplementaryComparativeModelRow is SupplementaryComparativeModelRow
    )
    assert (
        reports_api.SupplementaryComparativeModelTableResult
        is SupplementaryComparativeModelTableResult
    )
    assert (
        reports_api.write_supplementary_comparative_model_table
        is write_supplementary_comparative_model_table
    )
