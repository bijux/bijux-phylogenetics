from __future__ import annotations

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedNniMoveCandidate,
    iter_rooted_nni_move_candidates,
)


def test_topology_gateway_exports_rooted_nni_neighbor_surface() -> None:
    assert topology_api.RootedNniMoveCandidate is RootedNniMoveCandidate
    assert topology_api.iter_rooted_nni_move_candidates is iter_rooted_nni_move_candidates
