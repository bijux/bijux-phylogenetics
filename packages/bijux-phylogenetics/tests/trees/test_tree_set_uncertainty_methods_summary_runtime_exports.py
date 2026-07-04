from __future__ import annotations

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import write_tree_set_uncertainty_methods_summary_text


def test_tree_set_uncertainty_methods_summary_surfaces_export_publicly() -> None:
    assert callable(write_tree_set_uncertainty_methods_summary_text)
    assert (
        trees_api.write_tree_set_uncertainty_methods_summary_text
        is write_tree_set_uncertainty_methods_summary_text
    )
