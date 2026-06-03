from __future__ import annotations

import bijux_phylogenetics.comparative as comparative_api
from bijux_phylogenetics.comparative import write_comparative_methods_summary_text


def test_comparative_methods_summary_surfaces_export_publicly() -> None:
    assert callable(write_comparative_methods_summary_text)
    assert (
        comparative_api.write_comparative_methods_summary_text
        is write_comparative_methods_summary_text
    )
