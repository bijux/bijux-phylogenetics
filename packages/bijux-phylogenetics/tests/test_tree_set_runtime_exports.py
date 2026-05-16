from __future__ import annotations

import bijux_phylogenetics
from bijux_phylogenetics.tree_set import (
    TreeDistanceDistributionRow,
    TreeSetProcessingSummary,
    summarize_posterior_topology_diversity,
    write_tree_distance_distribution_table,
)


def test_public_runtime_exports_tree_set_scaling_surface() -> None:
    assert (
        bijux_phylogenetics.summarize_posterior_topology_diversity
        is summarize_posterior_topology_diversity
    )
    assert (
        bijux_phylogenetics.write_tree_distance_distribution_table
        is write_tree_distance_distribution_table
    )
    assert (
        bijux_phylogenetics.TreeDistanceDistributionRow is TreeDistanceDistributionRow
    )
    assert bijux_phylogenetics.TreeSetProcessingSummary is TreeSetProcessingSummary
