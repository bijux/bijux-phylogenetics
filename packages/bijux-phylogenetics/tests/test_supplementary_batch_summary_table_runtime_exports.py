from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
from bijux_phylogenetics.reports import (
    SupplementaryBatchSummaryRow,
    SupplementaryBatchSummaryTableResult,
    write_supplementary_batch_summary_table,
)


def test_supplementary_batch_summary_table_surfaces_export_publicly() -> None:
    assert reports_api.SupplementaryBatchSummaryRow is SupplementaryBatchSummaryRow
    assert (
        reports_api.SupplementaryBatchSummaryTableResult
        is SupplementaryBatchSummaryTableResult
    )
    assert (
        reports_api.write_supplementary_batch_summary_table
        is write_supplementary_batch_summary_table
    )
