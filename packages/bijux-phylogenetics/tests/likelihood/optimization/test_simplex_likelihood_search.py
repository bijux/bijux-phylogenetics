from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.phylo.likelihood.simplex_search import (
    resolve_simplex_coordinate_search_bounds,
    run_simplex_likelihood_search,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_simplex_likelihood_search_lifts_zero_starting_mass_and_keeps_candidates_valid() -> (
    None
):
    evaluations: list[dict[str, float]] = []

    def evaluate(candidate: dict[str, float]) -> tuple[dict[str, float], float]:
        evaluations.append(candidate)
        return candidate, (4.0 * candidate["B"]) - candidate["A"]

    result = run_simplex_likelihood_search(
        initial_values_by_name={"A": 0.999999, "B": 0.0, "C": 0.000001},
        expected_component_names=("A", "B", "C"),
        owner_name="simplex likelihood search test",
        evaluate=evaluate,
        minimum_component_value=1e-4,
        improvement_tolerance=1e-10,
        max_coordinate_passes=6,
    )

    assert result.initial_parameterization.constrained_mapping()["B"] >= 1e-4
    assert result.optimized_parameterization.constrained_mapping()["B"] > 0.9
    assert result.objective_value > 2.5
    assert result.function_evaluation_count > 1
    assert evaluations
    for candidate in evaluations:
        assert math.isclose(sum(candidate.values()), 1.0, rel_tol=0.0, abs_tol=1e-12)
        assert min(candidate.values()) >= 1e-4


def test_simplex_coordinate_search_bounds_follow_declared_floor() -> None:
    lower_bound, upper_bound = resolve_simplex_coordinate_search_bounds(
        4,
        minimum_component_value=1e-4,
        owner_name="simplex bound test",
    )

    assert math.isclose(lower_bound, -upper_bound, rel_tol=0.0, abs_tol=1e-12)
    assert upper_bound > 0.0


def test_simplex_coordinate_search_bounds_reject_floor_at_uniform_mass() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="minimum component value must stay below the uniform simplex mass",
    ):
        resolve_simplex_coordinate_search_bounds(
            4,
            minimum_component_value=0.25,
            owner_name="simplex bound test",
        )
