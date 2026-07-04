from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_invariant_proportion_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_invariant_proportion_move as propose_invariant_proportion_move_impl,
)


def test_bayesian_exports_invariant_proportion_proposal_surface() -> None:
    assert propose_invariant_proportion_move is propose_invariant_proportion_move_impl
