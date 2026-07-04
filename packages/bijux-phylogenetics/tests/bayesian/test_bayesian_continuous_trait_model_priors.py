from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_beta_continuous_trait_probability_prior,
    build_continuous_trait_model_prior_bundle,
    build_exponential_continuous_trait_scalar_prior,
    build_fixed_continuous_trait_scalar_prior,
    evaluate_continuous_trait_model_log_prior,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    ComparativeResidualSummary,
)
from bijux_phylogenetics.comparative.evolutionary_modes import (
    ContinuousEvolutionaryModeFitReport,
    fit_continuous_evolutionary_mode,
)
from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_geiger_continuous_fixture,
)


def _gamma_log_density(value: float, *, shape: float, scale: float) -> float:
    return (
        ((shape - 1.0) * math.log(value))
        - (value / scale)
        - math.lgamma(shape)
        - (shape * math.log(scale))
    )


def _lognormal_log_density(
    value: float,
    *,
    log_mean: float,
    log_standard_deviation: float,
) -> float:
    return (
        -math.log(value)
        - math.log(log_standard_deviation)
        - (0.5 * math.log(2.0 * math.pi))
        - (((math.log(value) - log_mean) ** 2) / (2.0 * (log_standard_deviation**2)))
    )


def _beta_log_density(value: float, *, alpha: float, beta: float) -> float:
    return (
        math.lgamma(alpha + beta)
        - math.lgamma(alpha)
        - math.lgamma(beta)
        + ((alpha - 1.0) * math.log(value))
        + ((beta - 1.0) * math.log1p(-value))
    )


def _build_synthetic_fit_report(
    *,
    mode: str,
    parameter_name: str | None,
    parameter_value: float | None,
    rate: float = 0.8,
    log_likelihood: float = -10.0,
) -> ContinuousEvolutionaryModeFitReport:
    return ContinuousEvolutionaryModeFitReport(
        tree_path=Path("synthetic-tree.nwk"),
        traits_path=Path("synthetic-traits.tsv"),
        taxon_column="taxon",
        trait="trait",
        taxon_count=4,
        taxa=["A", "B", "C", "D"],
        mode=mode,
        parameter_name=parameter_name,
        parameter_value=parameter_value,
        root_state=0.0,
        rate=rate,
        log_likelihood=log_likelihood,
        aic=0.0,
        aicc=0.0,
        likelihood_constant_policy="synthetic-policy",
        likelihood_comparison_policy="synthetic-comparison-policy",
        fitted_values=[0.0, 0.0, 0.0, 0.0],
        residuals=[0.0, 0.0, 0.0, 0.0],
        transformed_tree_newick="(A:1,B:1,(C:1,D:1):1);",
        confidence_intervals=[],
        residual_diagnostics=ComparativeResidualSummary(
            residual_mean=0.0,
            residual_variance=0.0,
            residual_skewness=0.0,
            max_abs_standardized_residual=0.0,
            phylogenetic_residual_lambda=0.0,
            outlier_taxa=[],
            warnings=[],
        ),
        optimizer_diagnostics=None,
        optimizer_profile_rows=None,
        identifiability_warnings=[],
        assumptions=[],
    )


def test_continuous_trait_model_prior_report_matches_analytical_ou_fixture() -> None:
    fit_report = _build_synthetic_fit_report(
        mode="ornstein-uhlenbeck",
        parameter_name="alpha",
        parameter_value=0.5,
        rate=0.8,
        log_likelihood=-10.0,
    )
    prior_bundle = build_continuous_trait_model_prior_bundle(
        alpha_prior=build_fixed_continuous_trait_scalar_prior(fixed_value=0.5),
        sigma_prior=build_exponential_continuous_trait_scalar_prior(rate=0.75),
    )

    report = evaluate_continuous_trait_model_log_prior(
        fit_report=fit_report,
        prior_bundle=prior_bundle,
    )

    expected_total_log_prior = 0.0 + (math.log(0.75) - (0.75 * 0.8))

    assert report.mode == "ornstein-uhlenbeck"
    assert report.prior_count == 2
    assert [row.target_name for row in report.rows] == ["alpha", "sigma"]
    assert report.log_likelihood == pytest.approx(-10.0)
    assert report.total_log_prior == pytest.approx(expected_total_log_prior, abs=1e-12)
    assert report.posterior_score == pytest.approx(
        report.log_likelihood + report.total_log_prior,
        abs=1e-12,
    )


@pytest.mark.parametrize(
    ("mode", "parameter_name", "parameter_value", "prior_bundle", "expected_targets"),
    [
        (
            "brownian",
            None,
            None,
            build_continuous_trait_model_prior_bundle(
                sigma_prior=build_exponential_continuous_trait_scalar_prior(rate=0.5)
            ),
            ["sigma"],
        ),
        (
            "pagel-lambda",
            "lambda",
            0.4,
            build_continuous_trait_model_prior_bundle(
                lambda_prior=build_beta_continuous_trait_probability_prior(
                    alpha=2.0,
                    beta=3.0,
                ),
                sigma_prior=build_exponential_continuous_trait_scalar_prior(rate=0.25),
            ),
            ["lambda", "sigma"],
        ),
        (
            "pagel-kappa",
            "kappa",
            0.0,
            build_continuous_trait_model_prior_bundle(
                kappa_prior=build_exponential_continuous_trait_scalar_prior(rate=1.1)
            ),
            ["kappa"],
        ),
        (
            "pagel-delta",
            "delta",
            1.25,
            build_continuous_trait_model_prior_bundle(
                delta_prior=build_fixed_continuous_trait_scalar_prior(fixed_value=1.25)
            ),
            ["delta"],
        ),
        (
            "ornstein-uhlenbeck",
            "alpha",
            0.0,
            build_continuous_trait_model_prior_bundle(
                alpha_prior=build_exponential_continuous_trait_scalar_prior(rate=0.9)
            ),
            ["alpha"],
        ),
        (
            "early-burst",
            "rate_change",
            0.0,
            build_continuous_trait_model_prior_bundle(
                rate_change_prior=build_exponential_continuous_trait_scalar_prior(
                    rate=0.6
                )
            ),
            ["rate_change"],
        ),
    ],
)
def test_continuous_trait_model_priors_support_owned_mode_targets(
    mode: str,
    parameter_name: str | None,
    parameter_value: float | None,
    prior_bundle: object,
    expected_targets: list[str],
) -> None:
    fit_report = _build_synthetic_fit_report(
        mode=mode,
        parameter_name=parameter_name,
        parameter_value=parameter_value,
    )

    report = evaluate_continuous_trait_model_log_prior(
        fit_report=fit_report,
        prior_bundle=prior_bundle,
    )

    assert [row.target_name for row in report.rows] == expected_targets
    assert report.posterior_score == pytest.approx(
        report.log_likelihood + report.total_log_prior,
        abs=1e-12,
    )


def test_continuous_trait_model_prior_changes_posterior_score_on_owned_fit_surface() -> (
    None
):
    fixture_entry = get_shared_geiger_continuous_fixture(
        "geiger_continuous_ou_known_truth_twenty_four_taxa"
    )
    fit_report = fit_continuous_evolutionary_mode(
        fixture_entry.tree_path,
        fixture_entry.traits_path,
        trait=fixture_entry.trait_name,
        taxon_column=fixture_entry.taxon_column,
        mode="ornstein-uhlenbeck",
    )
    first_prior_bundle = build_continuous_trait_model_prior_bundle(
        alpha_prior=build_exponential_continuous_trait_scalar_prior(rate=0.25),
        sigma_prior=build_exponential_continuous_trait_scalar_prior(rate=0.5),
    )
    second_prior_bundle = build_continuous_trait_model_prior_bundle(
        alpha_prior=build_exponential_continuous_trait_scalar_prior(rate=1.25),
        sigma_prior=build_exponential_continuous_trait_scalar_prior(rate=1.5),
    )

    first_report = evaluate_continuous_trait_model_log_prior(
        fit_report=fit_report,
        prior_bundle=first_prior_bundle,
    )
    second_report = evaluate_continuous_trait_model_log_prior(
        fit_report=fit_report,
        prior_bundle=second_prior_bundle,
    )

    assert first_report.log_likelihood == pytest.approx(fit_report.log_likelihood)
    assert second_report.log_likelihood == pytest.approx(fit_report.log_likelihood)
    assert not math.isclose(
        first_report.total_log_prior,
        second_report.total_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        first_report.posterior_score,
        second_report.posterior_score,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_continuous_trait_model_prior_beta_and_exponential_density_match_fixture() -> (
    None
):
    fit_report = _build_synthetic_fit_report(
        mode="pagel-lambda",
        parameter_name="lambda",
        parameter_value=0.4,
        rate=0.7,
        log_likelihood=-4.0,
    )
    prior_bundle = build_continuous_trait_model_prior_bundle(
        lambda_prior=build_beta_continuous_trait_probability_prior(
            alpha=2.5,
            beta=3.5,
        ),
        sigma_prior=build_exponential_continuous_trait_scalar_prior(rate=0.8),
    )

    report = evaluate_continuous_trait_model_log_prior(
        fit_report=fit_report,
        prior_bundle=prior_bundle,
    )

    expected_total_log_prior = _beta_log_density(0.4, alpha=2.5, beta=3.5) + (
        math.log(0.8) - (0.8 * 0.7)
    )

    assert report.total_log_prior == pytest.approx(expected_total_log_prior, abs=1e-12)
