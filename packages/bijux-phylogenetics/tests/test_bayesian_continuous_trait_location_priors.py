from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    CONTINUOUS_TRAIT_LOCATION_PRIOR_FAMILIES,
    build_fixed_continuous_trait_location_prior,
    build_normal_continuous_trait_location_prior,
    evaluate_continuous_trait_location_log_prior,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_normal_continuous_trait_location_prior_matches_analytical_density() -> None:
    prior_model = build_normal_continuous_trait_location_prior(
        mean=1.5,
        standard_deviation=0.75,
    )

    log_prior = evaluate_continuous_trait_location_log_prior(
        parameter_value=0.5,
        prior_model=prior_model,
        parameter_name="root_state",
    )

    centered_value = (0.5 - 1.5) / 0.75
    expected_log_prior = (
        -math.log(0.75) - (0.5 * math.log(2.0 * math.pi)) - (0.5 * (centered_value**2))
    )

    assert CONTINUOUS_TRAIT_LOCATION_PRIOR_FAMILIES == ("fixed", "normal")
    assert log_prior == pytest.approx(expected_log_prior, abs=1e-12)


def test_fixed_continuous_trait_location_prior_accepts_exact_value() -> None:
    prior_model = build_fixed_continuous_trait_location_prior(
        fixed_value=-0.25,
        fixed_tolerance=1e-9,
    )

    log_prior = evaluate_continuous_trait_location_log_prior(
        parameter_value=-0.25,
        prior_model=prior_model,
        parameter_name="root_state",
    )

    assert log_prior == 0.0


def test_fixed_continuous_trait_location_prior_rejects_out_of_tolerance_value() -> None:
    prior_model = build_fixed_continuous_trait_location_prior(
        fixed_value=0.0,
        fixed_tolerance=1e-6,
    )

    log_prior = evaluate_continuous_trait_location_log_prior(
        parameter_value=0.01,
        prior_model=prior_model,
        parameter_name="root_state",
    )

    assert log_prior == -math.inf


def test_normal_continuous_trait_location_prior_requires_positive_standard_deviation() -> (
    None
):
    with pytest.raises(
        PhylogeneticsError,
        match="requires standard_deviation > 0",
    ):
        build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=0.0,
        )
