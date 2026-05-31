from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_beta_probability_substitution_parameter_prior,
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_fixed_positive_substitution_parameter_prior,
    build_fixed_probability_substitution_parameter_prior,
    build_fixed_simplex_substitution_parameter_prior,
    build_lognormal_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
    evaluate_substitution_parameter_log_prior,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def _lognormal_log_density(
    value: float,
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> float:
    return (
        -math.log(value)
        - math.log(log_standard_deviation)
        - 0.5 * math.log(2.0 * math.pi)
        - ((math.log(value) - log_mean) ** 2) / (2.0 * (log_standard_deviation**2))
    )


def _beta_log_density(value: float, *, alpha: float, beta: float) -> float:
    return (
        math.lgamma(alpha + beta)
        - math.lgamma(alpha)
        - math.lgamma(beta)
        + ((alpha - 1.0) * math.log(value))
        + ((beta - 1.0) * math.log1p(-value))
    )


def _dirichlet_log_density(
    values: tuple[float, ...],
    concentration_parameters: tuple[float, ...],
) -> float:
    return (
        math.lgamma(math.fsum(concentration_parameters))
        - math.fsum(math.lgamma(alpha) for alpha in concentration_parameters)
        + math.fsum(
            (alpha - 1.0) * math.log(value)
            for alpha, value in zip(concentration_parameters, values, strict=True)
        )
    )


def test_substitution_parameter_prior_bundle_matches_analytical_density_fixture() -> (
    None
):
    prior_bundle = build_substitution_parameter_prior_bundle(
        kappa_prior=build_exponential_positive_substitution_parameter_prior(rate=0.75),
        exchangeability_prior=build_dirichlet_simplex_substitution_parameter_prior(
            expected_component_names=("AC", "AG", "AT", "CG", "CT", "GT"),
            concentration_parameters={
                "AC": 2.0,
                "AG": 3.0,
                "AT": 4.0,
                "CG": 5.0,
                "CT": 6.0,
                "GT": 7.0,
            },
        ),
        base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
            expected_component_names=("A", "C", "G", "T"),
            concentration_parameters={"A": 4.0, "C": 3.0, "G": 2.0, "T": 5.0},
        ),
        gamma_alpha_prior=build_lognormal_positive_substitution_parameter_prior(
            log_mean=math.log(0.8),
            log_standard_deviation=0.4,
        ),
        invariant_proportion_prior=build_beta_probability_substitution_parameter_prior(
            alpha=2.0,
            beta=5.0,
        ),
    )

    report = evaluate_substitution_parameter_log_prior(
        prior_bundle=prior_bundle,
        kappa=2.0,
        exchangeabilities={
            "AC": 1.0,
            "AG": 2.0,
            "AT": 1.0,
            "CG": 3.0,
            "CT": 2.0,
            "GT": 1.0,
        },
        base_frequencies={"A": 0.3, "C": 0.2, "G": 0.1, "T": 0.4},
        gamma_alpha=0.8,
        invariant_proportion=0.25,
    )

    expected_total_log_prior = (
        math.log(0.75)
        - (0.75 * 2.0)
        + _dirichlet_log_density(
            (0.1, 0.2, 0.1, 0.3, 0.2, 0.1),
            (2.0, 3.0, 4.0, 5.0, 6.0, 7.0),
        )
        + _dirichlet_log_density(
            (0.3, 0.2, 0.1, 0.4),
            (4.0, 3.0, 2.0, 5.0),
        )
        + _lognormal_log_density(
            0.8,
            log_mean=math.log(0.8),
            log_standard_deviation=0.4,
        )
        + _beta_log_density(0.25, alpha=2.0, beta=5.0)
    )

    assert report.prior_count == 5
    assert [row.target_name for row in report.rows] == [
        "kappa",
        "exchangeabilities",
        "base-frequencies",
        "gamma-alpha",
        "invariant-proportion",
    ]
    assert math.isclose(
        report.total_log_prior,
        expected_total_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_substitution_parameter_prior_bundle_changes_when_parameter_values_change() -> (
    None
):
    prior_bundle = build_substitution_parameter_prior_bundle(
        kappa_prior=build_exponential_positive_substitution_parameter_prior(rate=0.5),
        exchangeability_prior=build_dirichlet_simplex_substitution_parameter_prior(
            expected_component_names=("AC", "AG", "AT", "CG", "CT", "GT"),
            concentration_parameters=(2.0, 2.0, 2.0, 2.0, 2.0, 2.0),
        ),
        base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
            expected_component_names=("A", "C", "G", "T"),
            concentration_parameters=(3.0, 2.0, 2.0, 3.0),
        ),
    )

    first_report = evaluate_substitution_parameter_log_prior(
        prior_bundle=prior_bundle,
        kappa=1.5,
        exchangeabilities=(1.0, 2.0, 1.0, 3.0, 2.0, 1.0),
        base_frequencies=(0.3, 0.2, 0.1, 0.4),
    )
    second_report = evaluate_substitution_parameter_log_prior(
        prior_bundle=prior_bundle,
        kappa=2.5,
        exchangeabilities=(1.0, 1.0, 3.0, 1.0, 3.0, 1.0),
        base_frequencies=(0.2, 0.3, 0.2, 0.3),
    )

    assert not math.isclose(
        first_report.total_log_prior,
        second_report.total_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_substitution_parameter_fixed_priors_distinguish_exact_and_mismatched_values() -> (
    None
):
    prior_bundle = build_substitution_parameter_prior_bundle(
        kappa_prior=build_fixed_positive_substitution_parameter_prior(fixed_value=2.0),
        base_frequency_prior=build_fixed_simplex_substitution_parameter_prior(
            expected_component_names=("A", "C", "G", "T"),
            fixed_values={"A": 0.3, "C": 0.2, "G": 0.1, "T": 0.4},
        ),
        invariant_proportion_prior=build_fixed_probability_substitution_parameter_prior(
            fixed_value=0.25
        ),
    )

    exact_report = evaluate_substitution_parameter_log_prior(
        prior_bundle=prior_bundle,
        kappa=2.0,
        base_frequencies={"A": 0.3, "C": 0.2, "G": 0.1, "T": 0.4},
        invariant_proportion=0.25,
    )
    mismatched_report = evaluate_substitution_parameter_log_prior(
        prior_bundle=prior_bundle,
        kappa=2.2,
        base_frequencies={"A": 0.3, "C": 0.2, "G": 0.1, "T": 0.4},
        invariant_proportion=0.25,
    )

    assert exact_report.total_log_prior == pytest.approx(0.0, abs=1e-12)
    assert mismatched_report.total_log_prior == -math.inf


def test_substitution_parameter_prior_bundle_requires_values_for_configured_priors() -> (
    None
):
    prior_bundle = build_substitution_parameter_prior_bundle(
        gamma_alpha_prior=build_lognormal_positive_substitution_parameter_prior(
            log_mean=0.0,
            log_standard_deviation=0.5,
        )
    )

    with pytest.raises(
        PhylogeneticsError,
        match="requires gamma_alpha",
    ):
        evaluate_substitution_parameter_log_prior(prior_bundle=prior_bundle)
