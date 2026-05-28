from __future__ import annotations

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedTbrNeighborRow,
    RootedTbrNeighborhoodReport,
)


def test_topology_gateway_exports_rooted_tbr_neighbor_contracts() -> None:
    assert topology_api.RootedTbrNeighborRow is RootedTbrNeighborRow
    assert topology_api.RootedTbrNeighborhoodReport is RootedTbrNeighborhoodReport
