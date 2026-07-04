from __future__ import annotations

import bijux_phylogenetics.ancestral as ancestral_api
from bijux_phylogenetics.ancestral import write_ancestral_methods_summary_text


def test_ancestral_methods_summary_surfaces_export_publicly() -> None:
    assert callable(write_ancestral_methods_summary_text)
    assert (
        ancestral_api.write_ancestral_methods_summary_text
        is write_ancestral_methods_summary_text
    )
