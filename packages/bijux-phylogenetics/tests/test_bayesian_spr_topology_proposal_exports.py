from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_spr_topology_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_spr_topology_move as propose_spr_topology_move_impl,
)


def test_bayesian_exports_spr_topology_proposal_surface() -> None:
    assert propose_spr_topology_move is propose_spr_topology_move_impl
