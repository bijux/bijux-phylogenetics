from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_partition_linking_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_partition_linking_move as propose_partition_linking_move_impl,
)


def test_bayesian_exports_partition_linking_proposal_surface() -> None:
    assert propose_partition_linking_move is propose_partition_linking_move_impl
