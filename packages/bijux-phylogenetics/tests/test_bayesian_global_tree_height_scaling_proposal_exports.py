from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_global_tree_height_scaling_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_global_tree_height_scaling_move as propose_global_tree_height_scaling_move_impl,
)


def test_bayesian_exports_global_tree_height_scaling_proposal_surface() -> None:
    assert (
        propose_global_tree_height_scaling_move
        is propose_global_tree_height_scaling_move_impl
    )
