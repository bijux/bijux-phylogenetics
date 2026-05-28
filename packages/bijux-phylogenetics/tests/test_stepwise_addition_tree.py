from __future__ import annotations

import pytest

import bijux_phylogenetics.phylo.topology as topology_api
from bijux_phylogenetics.phylo.topology import (
    STEPWISE_ADDITION_ROOT_BRANCH_ID,
    validate_stepwise_addition_taxa,
    validate_stepwise_objective_direction,
)


def test_topology_gateway_exports_stepwise_addition_validation_surface() -> None:
    assert (
        topology_api.validate_stepwise_addition_taxa
        is validate_stepwise_addition_taxa
    )
    assert (
        topology_api.validate_stepwise_objective_direction
        is validate_stepwise_objective_direction
    )
    assert topology_api.STEPWISE_ADDITION_ROOT_BRANCH_ID == STEPWISE_ADDITION_ROOT_BRANCH_ID


def test_validate_stepwise_addition_taxa_preserves_insertion_order() -> None:
    assert validate_stepwise_addition_taxa(["Beta", "Alpha", "Gamma"]) == [
        "Beta",
        "Alpha",
        "Gamma",
    ]


@pytest.mark.parametrize(
    ("taxa", "message"),
    [
        (["Alpha"], "stepwise addition requires at least two taxa"),
        (
            ["Alpha", "", "Gamma"],
            "stepwise addition does not allow blank taxon labels",
        ),
        (
            ["Alpha", "Beta", "Alpha"],
            "stepwise addition requires distinct taxa; duplicates: Alpha",
        ),
    ],
)
def test_validate_stepwise_addition_taxa_rejects_invalid_taxon_sets(
    taxa: list[str],
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        validate_stepwise_addition_taxa(taxa)


def test_validate_stepwise_objective_direction_accepts_supported_values() -> None:
    assert validate_stepwise_objective_direction(" minimize ") == "minimize"
    assert validate_stepwise_objective_direction("MAXIMIZE") == "maximize"


def test_validate_stepwise_objective_direction_rejects_unknown_values() -> None:
    with pytest.raises(
        ValueError,
        match="objective_direction must be one of maximize, minimize",
    ):
        validate_stepwise_objective_direction("median")
