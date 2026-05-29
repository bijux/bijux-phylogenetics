from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_discrete_trait_rate_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_discrete_trait_rate_move as propose_discrete_trait_rate_move_impl,
)


def test_bayesian_exports_discrete_trait_rate_proposal_surface() -> None:
    assert propose_discrete_trait_rate_move is propose_discrete_trait_rate_move_impl
