from __future__ import annotations

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    TopologyNeighborhoodSummaryReport,
    summarize_topology_neighborhood,
    write_topology_neighborhood_summary_table,
)


def test_public_runtime_exports_topology_neighborhood_summary_surface() -> None:
    assert (
        topology_api.TopologyNeighborhoodSummaryReport
        is TopologyNeighborhoodSummaryReport
    )
    assert (
        topology_api.summarize_topology_neighborhood is summarize_topology_neighborhood
    )
    assert (
        topology_api.write_topology_neighborhood_summary_table
        is write_topology_neighborhood_summary_table
    )
