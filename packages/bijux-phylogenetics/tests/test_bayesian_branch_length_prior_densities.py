from __future__ import annotations

import math

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
    build_fixed_branch_length_prior,
    build_gamma_branch_length_prior,
    build_lognormal_branch_length_prior,
    evaluate_branch_length_log_prior,
)


def test_exponential_branch_length_prior_matches_analytical_fixture() -> None:
    prior_model = build_exponential_branch_length_prior(rate=5.0)

    observed = evaluate_branch_length_log_prior(0.2, prior_model)
    expected = math.log(5.0) - (5.0 * 0.2)

    assert math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-12)


def test_gamma_branch_length_prior_matches_analytical_fixture() -> None:
    prior_model = build_gamma_branch_length_prior(shape=3.0, scale=0.5)

    observed = evaluate_branch_length_log_prior(1.25, prior_model)
    expected = (
        (3.0 - 1.0) * math.log(1.25)
        - (1.25 / 0.5)
        - math.lgamma(3.0)
        - (3.0 * math.log(0.5))
    )

    assert math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-12)


def test_lognormal_branch_length_prior_matches_analytical_fixture() -> None:
    prior_model = build_lognormal_branch_length_prior(
        log_mean=-1.0,
        log_standard_deviation=0.5,
    )

    observed = evaluate_branch_length_log_prior(0.25, prior_model)
    centered = math.log(0.25) - (-1.0)
    expected = (
        -math.log(0.25)
        - math.log(0.5)
        - (0.5 * math.log(2.0 * math.pi))
        - (centered * centered / (2.0 * 0.25))
    )

    assert math.isclose(observed, expected, rel_tol=0.0, abs_tol=1e-12)


def test_fixed_branch_length_prior_matches_and_rejects_mismatch() -> None:
    prior_model = build_fixed_branch_length_prior(
        fixed_value=0.1,
        fixed_tolerance=1e-9,
    )

    assert evaluate_branch_length_log_prior(0.1, prior_model) == 0.0
    assert evaluate_branch_length_log_prior(0.1000000005, prior_model) == 0.0
    assert evaluate_branch_length_log_prior(0.11, prior_model) == -math.inf
