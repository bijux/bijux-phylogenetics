from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_base_frequency_simplex_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_base_frequency_simplex_move as propose_base_frequency_simplex_move_impl,
)


def test_bayesian_exports_base_frequency_simplex_proposal_surface() -> None:
    assert (
        propose_base_frequency_simplex_move is propose_base_frequency_simplex_move_impl
    )
