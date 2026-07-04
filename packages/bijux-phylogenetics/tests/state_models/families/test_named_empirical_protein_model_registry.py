from __future__ import annotations

import math

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood import (
    empirical_protein_exchangeability_matrix,
    empirical_protein_rate_matrix,
    empirical_protein_stationary_frequencies,
    list_supported_empirical_protein_model_names,
    validate_empirical_protein_model_name,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    compute_ctmc_expected_substitution_rate,
    solve_ctmc_stationary_distribution,
)
from bijux_phylogenetics.phylo.likelihood.protein import PROTEIN_STATE_ORDER


def test_supported_empirical_protein_model_names_are_explicit() -> None:
    assert list_supported_empirical_protein_model_names() == ("JTT", "WAG", "LG")
    assert validate_empirical_protein_model_name(" jtt ") == "JTT"
    assert validate_empirical_protein_model_name("wag") == "WAG"
    assert validate_empirical_protein_model_name("Lg") == "LG"


def test_empirical_protein_model_registry_rejects_unknown_name() -> None:
    with pytest.raises(ValueError, match="empirical protein model name must be one of"):
        validate_empirical_protein_model_name("q.insect")


@pytest.mark.parametrize("model_name", ("JTT", "WAG", "LG"))
def test_named_empirical_protein_models_stay_normalized_and_stationary(
    model_name: str,
) -> None:
    exchangeability_matrix = empirical_protein_exchangeability_matrix(model_name)
    stationary_frequencies = empirical_protein_stationary_frequencies(model_name)
    rate_matrix = empirical_protein_rate_matrix(model_name)
    stationary_report = solve_ctmc_stationary_distribution(
        rate_matrix,
        state_labels=PROTEIN_STATE_ORDER,
    )

    assert exchangeability_matrix.shape == (20, 20)
    assert rate_matrix.shape == (20, 20)
    assert stationary_frequencies.shape == (20,)
    assert numpy.allclose(
        exchangeability_matrix,
        exchangeability_matrix.T,
        rtol=0.0,
        atol=1e-12,
    )
    assert numpy.all(exchangeability_matrix >= 0.0)
    assert numpy.all(stationary_frequencies > 0.0)
    assert math.isclose(
        float(stationary_frequencies.sum()),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert numpy.all(rate_matrix[numpy.eye(20, dtype=bool)] < 0.0)
    assert numpy.all(rate_matrix[~numpy.eye(20, dtype=bool)] >= 0.0)
    assert numpy.allclose(
        rate_matrix.sum(axis=1),
        numpy.zeros(20, dtype=float),
        rtol=0.0,
        atol=1e-12,
    )
    assert math.isclose(
        compute_ctmc_expected_substitution_rate(
            rate_matrix,
            stationary_frequencies,
            state_labels=PROTEIN_STATE_ORDER,
        ),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert numpy.allclose(
        stationary_report.probabilities,
        stationary_frequencies,
        rtol=0.0,
        atol=1e-8,
    )
