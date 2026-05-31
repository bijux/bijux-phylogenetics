from __future__ import annotations

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    TREE_SET_SPLIT_FREQUENCY_POLICIES,
    TreeSetSplitFrequencyReport,
    TreeSetSplitFrequencyRow,
    compute_tree_set_split_frequency_table,
    write_tree_set_split_frequency_table,
)


def test_public_runtime_exports_tree_set_split_frequency_surface() -> None:
    assert (
        trees_api.TREE_SET_SPLIT_FREQUENCY_POLICIES is TREE_SET_SPLIT_FREQUENCY_POLICIES
    )
    assert trees_api.TreeSetSplitFrequencyRow is TreeSetSplitFrequencyRow
    assert trees_api.TreeSetSplitFrequencyReport is TreeSetSplitFrequencyReport
    assert (
        trees_api.compute_tree_set_split_frequency_table
        is compute_tree_set_split_frequency_table
    )
    assert (
        trees_api.write_tree_set_split_frequency_table
        is write_tree_set_split_frequency_table
    )
