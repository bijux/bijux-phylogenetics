from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
from bijux_phylogenetics.reports import (
    SupplementaryModelSelectionRow,
    SupplementaryModelSelectionTableResult,
    write_supplementary_model_selection_table,
)


def test_supplementary_model_selection_table_surfaces_export_publicly() -> None:
    assert reports_api.SupplementaryModelSelectionRow is SupplementaryModelSelectionRow
    assert (
        reports_api.SupplementaryModelSelectionTableResult
        is SupplementaryModelSelectionTableResult
    )
    assert (
        reports_api.write_supplementary_model_selection_table
        is write_supplementary_model_selection_table
    )
