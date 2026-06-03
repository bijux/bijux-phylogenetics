from __future__ import annotations

import pytest

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    RootedSprEnumerationBudget,
    validate_rooted_spr_enumeration_budget,
)


def test_topology_gateway_exports_rooted_spr_budget_surface() -> None:
    assert topology_api.RootedSprEnumerationBudget is RootedSprEnumerationBudget
    assert (
        topology_api.validate_rooted_spr_enumeration_budget
        is validate_rooted_spr_enumeration_budget
    )


def test_rooted_spr_budget_validation_accepts_empty_budget() -> None:
    budget = validate_rooted_spr_enumeration_budget(None)

    assert budget == RootedSprEnumerationBudget()


def test_rooted_spr_budget_validation_accepts_positive_limits() -> None:
    budget = validate_rooted_spr_enumeration_budget(
        RootedSprEnumerationBudget(
            max_pruned_clade_count=2,
            max_regraft_target_count_per_pruned_clade=5,
        )
    )

    assert budget.max_pruned_clade_count == 2
    assert budget.max_regraft_target_count_per_pruned_clade == 5


def test_rooted_spr_budget_validation_rejects_nonpositive_prune_limit() -> None:
    with pytest.raises(
        ValueError,
        match="rooted SPR prune-node budget must be positive when provided",
    ):
        validate_rooted_spr_enumeration_budget(
            RootedSprEnumerationBudget(max_pruned_clade_count=0)
        )


def test_rooted_spr_budget_validation_rejects_nonpositive_regraft_limit() -> None:
    with pytest.raises(
        ValueError,
        match="rooted SPR regraft-target budget must be positive when provided",
    ):
        validate_rooted_spr_enumeration_budget(
            RootedSprEnumerationBudget(max_regraft_target_count_per_pruned_clade=0)
        )
