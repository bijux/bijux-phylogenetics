from __future__ import annotations

from bijux_phylogenetics.reports import (
    TreeInferenceMethodsSummaryTextResult,
    write_tree_inference_methods_summary_text,
)


def test_tree_inference_methods_summary_runtime_exports() -> None:
    assert TreeInferenceMethodsSummaryTextResult.__name__ == (
        "TreeInferenceMethodsSummaryTextResult"
    )
    assert callable(write_tree_inference_methods_summary_text)
