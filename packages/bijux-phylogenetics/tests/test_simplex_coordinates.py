from __future__ import annotations

import math

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood.dna import (
    normalize_dna_exchangeabilities_by_anchor,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    DNA_EXCHANGEABILITY_LABELS,
    parameterize_dna_base_frequency_simplex,
    parameterize_dna_exchangeability_simplex,
    resolve_anchor_normalized_dna_exchangeabilities_from_unconstrained,
    resolve_dna_base_frequencies_from_unconstrained,
    resolve_dna_exchangeability_simplex_from_unconstrained,
)
from bijux_phylogenetics.phylo.likelihood.simplex_coordinates import (
    parameterize_named_simplex,
    resolve_named_simplex_from_unconstrained,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_named_simplex_round_trip_preserves_component_names_and_values() -> None:
    parameterization = parameterize_named_simplex(
        {"left": 0.2, "middle": 0.3, "right": 0.5},
        expected_component_names=("left", "middle", "right"),
        owner_name="generic simplex test",
        reference_component_name="right",
    )

    round_trip = resolve_named_simplex_from_unconstrained(
        parameterization.unconstrained_values,
        component_names=("left", "middle", "right"),
        owner_name="generic simplex test",
        reference_component_name="right",
    )

    assert round_trip.component_names == ("left", "middle", "right")
    assert round_trip.reference_component_name == "right"
    assert round_trip.constrained_mapping() == pytest.approx(
        {"left": 0.2, "middle": 0.3, "right": 0.5}
    )


def test_dna_base_frequency_simplex_round_trip_preserves_values() -> None:
    parameterization = parameterize_dna_base_frequency_simplex(
        {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
    )

    round_trip = resolve_dna_base_frequencies_from_unconstrained(
        parameterization.unconstrained_values
    )

    assert parameterization.component_names == ("A", "C", "G", "T")
    assert parameterization.reference_component_name == "T"
    assert parameterization.constrained_mapping() == pytest.approx(
        {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
    )
    assert numpy.allclose(
        round_trip,
        numpy.array([0.4, 0.1, 0.2, 0.3], dtype=float),
        rtol=0.0,
        atol=1e-12,
    )
    assert math.isclose(float(round_trip.sum()), 1.0, rel_tol=0.0, abs_tol=1e-12)


def test_dna_exchangeability_simplex_round_trip_preserves_normalized_values() -> None:
    anchored_exchangeabilities = {
        "AC": 1.0,
        "AG": 10.436093496916191,
        "AT": 0.624578992788826,
        "CG": 2.351009315272331,
        "CT": 2.066395062525222,
        "GT": 5.880958916954807,
    }
    parameterization = parameterize_dna_exchangeability_simplex(
        anchored_exchangeabilities
    )

    round_trip = resolve_dna_exchangeability_simplex_from_unconstrained(
        parameterization.unconstrained_values
    )

    assert parameterization.component_names == DNA_EXCHANGEABILITY_LABELS
    assert parameterization.reference_component_name == "GT"
    assert round_trip.constrained_values == pytest.approx(
        parameterization.constrained_values,
        rel=0.0,
        abs=1e-12,
    )
    assert math.isclose(
        math.fsum(round_trip.constrained_values),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_dna_exchangeability_simplex_round_trip_preserves_anchor_normalized_values() -> (
    None
):
    anchored_exchangeabilities = {
        "AC": 1.0,
        "AG": 10.436093496916191,
        "AT": 0.624578992788826,
        "CG": 2.351009315272331,
        "CT": 2.066395062525222,
        "GT": 5.880958916954807,
    }
    parameterization = parameterize_dna_exchangeability_simplex(
        anchored_exchangeabilities
    )

    round_trip = resolve_anchor_normalized_dna_exchangeabilities_from_unconstrained(
        parameterization.unconstrained_values
    )

    assert numpy.allclose(
        round_trip,
        normalize_dna_exchangeabilities_by_anchor(
            anchored_exchangeabilities,
            model_name="GTR",
        ),
        rtol=0.0,
        atol=1e-12,
    )


def test_dna_base_frequency_simplex_rejects_boundary_component() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        parameterize_dna_base_frequency_simplex(
            {"A": 0.5, "C": 0.5, "G": 0.0, "T": 0.0}
        )

    assert error_info.value.code == "simplex_coordinate_component_not_positive"
    assert error_info.value.details["component_name"] == "G"
    assert error_info.value.details["component_value"] == 0.0


def test_dna_base_frequency_simplex_rejects_wrong_unconstrained_dimension() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        resolve_dna_base_frequencies_from_unconstrained([0.1, 0.2])

    assert error_info.value.code == "simplex_coordinate_dimension_mismatch"
    assert error_info.value.details["expected_dimension"] == 3
    assert error_info.value.details["observed_dimension"] == 2
