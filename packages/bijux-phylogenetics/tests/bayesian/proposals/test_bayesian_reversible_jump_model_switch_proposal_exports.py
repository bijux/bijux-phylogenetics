from __future__ import annotations

from bijux_phylogenetics.bayesian import propose_reversible_jump_model_switch_move
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_reversible_jump_model_switch_move as propose_reversible_jump_model_switch_move_impl,
)


def test_bayesian_exports_reversible_jump_model_switch_proposal_surface() -> None:
    assert (
        propose_reversible_jump_model_switch_move
        is propose_reversible_jump_model_switch_move_impl
    )
