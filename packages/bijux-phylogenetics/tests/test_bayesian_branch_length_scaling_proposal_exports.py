from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_branch_length_scaling_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_branch_length_scaling_move as propose_branch_length_scaling_move_impl,
)


def test_bayesian_exports_branch_length_scaling_proposal_surface() -> None:
    assert propose_branch_length_scaling_move is propose_branch_length_scaling_move_impl
