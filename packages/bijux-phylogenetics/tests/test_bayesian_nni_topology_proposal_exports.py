from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_nni_topology_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_nni_topology_move as propose_nni_topology_move_impl,
)


def test_bayesian_exports_nni_topology_proposal_surface() -> None:
    assert propose_nni_topology_move is propose_nni_topology_move_impl
