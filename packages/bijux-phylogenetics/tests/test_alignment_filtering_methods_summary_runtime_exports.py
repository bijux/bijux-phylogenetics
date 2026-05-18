from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
from bijux_phylogenetics.reports import (
    AlignmentFilteringMethodsSummaryTextResult,
    write_alignment_filtering_methods_summary_text,
)


def test_alignment_filtering_methods_summary_surfaces_export_publicly() -> None:
    assert (
        reports_api.AlignmentFilteringMethodsSummaryTextResult
        is AlignmentFilteringMethodsSummaryTextResult
    )
    assert (
        reports_api.write_alignment_filtering_methods_summary_text
        is write_alignment_filtering_methods_summary_text
    )
