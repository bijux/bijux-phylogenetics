from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_gamma_alpha_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_gamma_alpha_move as propose_gamma_alpha_move_impl,
)


def test_bayesian_exports_gamma_alpha_proposal_surface() -> None:
    assert propose_gamma_alpha_move is propose_gamma_alpha_move_impl
