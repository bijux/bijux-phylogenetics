from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
from bijux_phylogenetics.reports import (
    SupplementaryAlignmentDiagnosticsRow,
    SupplementaryAlignmentDiagnosticsTableResult,
    write_supplementary_alignment_diagnostics_table,
)


def test_supplementary_alignment_table_surfaces_export_publicly() -> None:
    assert (
        reports_api.SupplementaryAlignmentDiagnosticsRow
        is SupplementaryAlignmentDiagnosticsRow
    )
    assert (
        reports_api.SupplementaryAlignmentDiagnosticsTableResult
        is SupplementaryAlignmentDiagnosticsTableResult
    )
    assert (
        reports_api.write_supplementary_alignment_diagnostics_table
        is write_supplementary_alignment_diagnostics_table
    )
