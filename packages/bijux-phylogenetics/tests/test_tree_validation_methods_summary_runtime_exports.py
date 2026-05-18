from __future__ import annotations

import bijux_phylogenetics.reports as reports_api
from bijux_phylogenetics.reports import (
    TreeValidationMethodsSummaryTextResult,
    write_tree_validation_methods_summary_text,
)


def test_tree_validation_methods_summary_surfaces_export_publicly() -> None:
    assert (
        reports_api.TreeValidationMethodsSummaryTextResult
        is TreeValidationMethodsSummaryTextResult
    )
    assert (
        reports_api.write_tree_validation_methods_summary_text
        is write_tree_validation_methods_summary_text
    )
