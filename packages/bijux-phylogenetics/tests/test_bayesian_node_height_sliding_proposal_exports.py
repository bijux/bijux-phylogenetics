from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_node_height_sliding_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_node_height_sliding_move as propose_node_height_sliding_move_impl,
)


def test_bayesian_exports_node_height_sliding_proposal_surface() -> None:
    assert propose_node_height_sliding_move is propose_node_height_sliding_move_impl
