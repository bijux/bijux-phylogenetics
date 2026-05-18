from __future__ import annotations

import bijux_phylogenetics.comparative as comparative_api
from bijux_phylogenetics.comparative import write_diversification_methods_summary_text


def test_diversification_methods_summary_surfaces_export_publicly() -> None:
    assert callable(write_diversification_methods_summary_text)
    assert (
        comparative_api.write_diversification_methods_summary_text
        is write_diversification_methods_summary_text
    )
