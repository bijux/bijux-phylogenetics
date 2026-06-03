from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
from bijux_phylogenetics.reports import (
    SupplementaryTreeDiagnosticsRow,
    SupplementaryTreeDiagnosticsTableResult,
    write_supplementary_tree_diagnostics_table,
)


def test_supplementary_tree_table_surfaces_export_publicly() -> None:
    assert (
        reports_api.SupplementaryTreeDiagnosticsRow is SupplementaryTreeDiagnosticsRow
    )
    assert (
        reports_api.SupplementaryTreeDiagnosticsTableResult
        is SupplementaryTreeDiagnosticsTableResult
    )
    assert (
        reports_api.write_supplementary_tree_diagnostics_table
        is write_supplementary_tree_diagnostics_table
    )
