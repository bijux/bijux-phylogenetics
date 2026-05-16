from __future__ import annotations

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    TreeDistanceDistributionRow,
    TreeSetProcessingSummary,
    summarize_posterior_topology_diversity,
    write_tree_distance_distribution_table,
)


def test_public_runtime_exports_tree_set_scaling_surface() -> None:
    assert (
        trees_api.summarize_posterior_topology_diversity
        is summarize_posterior_topology_diversity
    )
    assert (
        trees_api.write_tree_distance_distribution_table
        is write_tree_distance_distribution_table
    )
    assert trees_api.TreeDistanceDistributionRow is TreeDistanceDistributionRow
    assert trees_api.TreeSetProcessingSummary is TreeSetProcessingSummary
