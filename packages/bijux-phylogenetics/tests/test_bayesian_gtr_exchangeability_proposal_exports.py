from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_gtr_exchangeability_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_gtr_exchangeability_move as propose_gtr_exchangeability_move_impl,
)


def test_bayesian_exports_gtr_exchangeability_proposal_surface() -> None:
    assert propose_gtr_exchangeability_move is propose_gtr_exchangeability_move_impl
