from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.bayesian.clock_model_priors import (
    CLOCK_MODEL_SCALAR_PRIOR_FAMILIES,
    ClockModelScalarPriorModel,
    build_exponential_clock_model_scalar_prior,
    build_fixed_clock_model_scalar_prior,
    build_gamma_clock_model_scalar_prior,
    build_lognormal_clock_model_scalar_prior,
    evaluate_clock_model_scalar_log_prior,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_clock_model_scalar_prior_builders_preserve_hyperparameters() -> None:
    exponential_prior = build_exponential_clock_model_scalar_prior(rate=2.0)
    gamma_prior = build_gamma_clock_model_scalar_prior(shape=3.0, scale=0.4)
    lognormal_prior = build_lognormal_clock_model_scalar_prior(
        log_mean=-1.0,
        log_standard_deviation=0.5,
    )
    fixed_prior = build_fixed_clock_model_scalar_prior(fixed_value=0.25)

    assert CLOCK_MODEL_SCALAR_PRIOR_FAMILIES == (
        "exponential",
        "fixed",
        "gamma",
        "lognormal",
    )
    assert exponential_prior.parameter_values() == {"rate": 2.0}
    assert gamma_prior.parameter_values() == {"shape": 3.0, "scale": 0.4}
    assert lognormal_prior.parameter_values() == {
        "log_mean": -1.0,
        "log_standard_deviation": 0.5,
    }
    assert fixed_prior.parameter_values() == {
        "fixed_value": 0.25,
        "fixed_tolerance": 1e-12,
    }


def test_clock_model_scalar_prior_evaluation_matches_expected_densities() -> None:
    exponential_log_prior = evaluate_clock_model_scalar_log_prior(
        parameter_value=0.4,
        prior_model=build_exponential_clock_model_scalar_prior(rate=2.0),
        parameter_name="mean-clock-rate",
    )
    gamma_log_prior = evaluate_clock_model_scalar_log_prior(
        parameter_value=0.4,
        prior_model=build_gamma_clock_model_scalar_prior(shape=3.0, scale=0.4),
        parameter_name="mean-clock-rate",
    )
    lognormal_log_prior = evaluate_clock_model_scalar_log_prior(
        parameter_value=0.4,
        prior_model=build_lognormal_clock_model_scalar_prior(
            log_mean=-1.0,
            log_standard_deviation=0.5,
        ),
        parameter_name="log-standard-deviation",
    )
    fixed_log_prior = evaluate_clock_model_scalar_log_prior(
        parameter_value=0.25,
        prior_model=build_fixed_clock_model_scalar_prior(fixed_value=0.25),
        parameter_name="mean-clock-rate",
    )

    assert math.isclose(
        exponential_log_prior,
        math.log(2.0) - (2.0 * 0.4),
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        gamma_log_prior,
        ((3.0 - 1.0) * math.log(0.4))
        - (0.4 / 0.4)
        - math.lgamma(3.0)
        - (3.0 * math.log(0.4)),
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        lognormal_log_prior,
        -math.log(0.4)
        - math.log(0.5)
        - (0.5 * math.log(2.0 * math.pi))
        - (((math.log(0.4) + 1.0) ** 2) / (2.0 * (0.5**2))),
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert fixed_log_prior == 0.0


def test_clock_model_fixed_prior_returns_negative_infinity_outside_tolerance() -> None:
    log_prior = evaluate_clock_model_scalar_log_prior(
        parameter_value=0.4,
        prior_model=build_fixed_clock_model_scalar_prior(
            fixed_value=0.25,
            fixed_tolerance=1e-6,
        ),
        parameter_name="mean-clock-rate",
    )

    assert log_prior == -math.inf


@pytest.mark.parametrize(
    ("builder", "kwargs", "message"),
    [
        (
            build_exponential_clock_model_scalar_prior,
            {"rate": 0.0},
            "requires 'rate' to be positive and finite",
        ),
        (
            build_gamma_clock_model_scalar_prior,
            {"shape": 0.0, "scale": 1.0},
            "requires 'shape' to be positive and finite",
        ),
        (
            build_gamma_clock_model_scalar_prior,
            {"shape": 2.0, "scale": 0.0},
            "requires 'scale' to be positive and finite",
        ),
        (
            build_lognormal_clock_model_scalar_prior,
            {"log_mean": 0.0, "log_standard_deviation": 0.0},
            "requires 'log_standard_deviation' to be positive and finite",
        ),
        (
            build_fixed_clock_model_scalar_prior,
            {"fixed_value": 0.0},
            "requires 'fixed_value' to be positive and finite",
        ),
    ],
)
def test_clock_model_scalar_prior_builders_reject_invalid_parameters(
    builder,
    kwargs,
    message: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=message):
        builder(**kwargs)


def test_clock_model_scalar_prior_evaluation_requires_positive_parameter_values() -> (
    None
):
    prior_model = build_exponential_clock_model_scalar_prior(rate=2.0)

    with pytest.raises(
        PhylogeneticsError,
        match="requires 'mean-clock-rate' to be positive and finite",
    ):
        evaluate_clock_model_scalar_log_prior(
            parameter_value=0.0,
            prior_model=prior_model,
            parameter_name="mean-clock-rate",
        )


def test_clock_model_scalar_prior_evaluation_rejects_unsupported_family() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="clock-model scalar prior family is unsupported",
    ):
        evaluate_clock_model_scalar_log_prior(
            parameter_value=0.4,
            prior_model=ClockModelScalarPriorModel(family="unsupported"),
            parameter_name="mean-clock-rate",
        )
