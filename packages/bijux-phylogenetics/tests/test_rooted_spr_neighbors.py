from __future__ import annotations

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedSprMoveCandidate,
    RootedSprNeighborRow,
    RootedSprNeighborhoodReport,
)


def test_topology_gateway_exports_rooted_spr_neighbor_contracts() -> None:
    assert topology_api.RootedSprMoveCandidate is RootedSprMoveCandidate
    assert topology_api.RootedSprNeighborRow is RootedSprNeighborRow
    assert topology_api.RootedSprNeighborhoodReport is RootedSprNeighborhoodReport
