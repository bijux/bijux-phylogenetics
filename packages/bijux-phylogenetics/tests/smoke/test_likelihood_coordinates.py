from __future__ import annotations

import math

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    parameterize_dna_base_frequency_simplex,
    resolve_dna_base_frequencies_from_unconstrained,
)
from bijux_phylogenetics.phylo.likelihood.simplex_coordinates import (
    parameterize_named_simplex,
    resolve_named_simplex_from_unconstrained,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

pytestmark = pytest.mark.smoke


def test_named_simplex_round_trip_preserves_component_mapping() -> None:
    parameterization = parameterize_named_simplex(
        {"left": 0.2, "middle": 0.3, "right": 0.5},
        expected_component_names=("left", "middle", "right"),
        owner_name="smoke simplex",
        reference_component_name="right",
    )

    round_trip = resolve_named_simplex_from_unconstrained(
        parameterization.unconstrained_values,
        component_names=("left", "middle", "right"),
        owner_name="smoke simplex",
        reference_component_name="right",
    )

    assert round_trip.constrained_mapping() == pytest.approx(
        {"left": 0.2, "middle": 0.3, "right": 0.5}
    )


def test_dna_base_frequency_simplex_round_trip_preserves_total_mass() -> None:
    parameterization = parameterize_dna_base_frequency_simplex(
        {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
    )

    round_trip = resolve_dna_base_frequencies_from_unconstrained(
        parameterization.unconstrained_values
    )

    assert numpy.allclose(
        round_trip,
        numpy.array([0.4, 0.1, 0.2, 0.3], dtype=float),
        rtol=0.0,
        atol=1e-12,
    )
    assert math.isclose(float(round_trip.sum()), 1.0, rel_tol=0.0, abs_tol=1e-12)


def test_dna_base_frequency_simplex_rejects_boundary_component() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        parameterize_dna_base_frequency_simplex(
            {"A": 0.5, "C": 0.5, "G": 0.0, "T": 0.0}
        )

    assert error_info.value.code == "simplex_coordinate_component_not_positive"
