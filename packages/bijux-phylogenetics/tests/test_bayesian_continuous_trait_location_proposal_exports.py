from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_continuous_trait_location_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_continuous_trait_location_move as propose_continuous_trait_location_move_impl,
)


def test_bayesian_exports_continuous_trait_location_proposal() -> None:
    assert (
        propose_continuous_trait_location_move
        is propose_continuous_trait_location_move_impl
    )
