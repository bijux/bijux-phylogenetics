from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.calibration_priors import (
    load_calibration_prior_definitions,
)
from bijux_phylogenetics.bayesian.clock_model_priors import (
    build_gamma_clock_model_scalar_prior,
)
from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    build_normal_continuous_trait_location_prior,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_beta_continuous_trait_probability_prior,
    build_lognormal_continuous_trait_scalar_prior,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    build_exponential_discrete_trait_rate_prior,
)
from bijux_phylogenetics.bayesian.prior_protocol import (
    BayesianPriorModelSummary,
    summarize_bayesian_prior_model,
    supports_bayesian_prior_model_protocol,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_beta_probability_substitution_parameter_prior,
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
)
from bijux_phylogenetics.bayesian.time_tree_priors import (
    SkylineCoalescentEpoch,
    build_constant_population_coalescent_tree_prior,
    build_crown_conditioned_birth_death_tree_prior,
    build_crown_conditioned_yule_tree_prior,
    build_skyline_coalescent_tree_prior,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


@pytest.mark.parametrize(
    ("prior_model", "expected_family", "expected_parameter_values"),
    [
        (
            build_exponential_branch_length_prior(rate=2.0),
            "exponential",
            {"rate": 2.0},
        ),
        (
            build_exponential_positive_substitution_parameter_prior(rate=0.75),
            "exponential",
            {"rate": 0.75},
        ),
        (
            build_beta_probability_substitution_parameter_prior(
                alpha=2.5,
                beta=3.5,
            ),
            "beta",
            {"alpha": 2.5, "beta": 3.5},
        ),
        (
            build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=("A", "C", "G", "T"),
                concentration_parameters=(1.0, 2.0, 3.0, 4.0),
            ),
            "dirichlet",
            {"alpha_A": 1.0, "alpha_C": 2.0, "alpha_G": 3.0, "alpha_T": 4.0},
        ),
        (
            build_gamma_clock_model_scalar_prior(shape=2.0, scale=0.5),
            "gamma",
            {"shape": 2.0, "scale": 0.5},
        ),
        (
            build_exponential_discrete_trait_rate_prior(rate=1.25),
            "exponential",
            {"rate": 1.25},
        ),
        (
            build_normal_continuous_trait_location_prior(
                mean=0.25,
                standard_deviation=1.5,
            ),
            "normal",
            {"mean": 0.25, "standard_deviation": 1.5},
        ),
        (
            build_lognormal_continuous_trait_scalar_prior(
                log_mean=-0.5,
                log_standard_deviation=0.75,
            ),
            "lognormal",
            {"log_mean": -0.5, "log_standard_deviation": 0.75},
        ),
        (
            build_beta_continuous_trait_probability_prior(alpha=1.5, beta=2.5),
            "beta",
            {"alpha": 1.5, "beta": 2.5},
        ),
        (
            build_crown_conditioned_yule_tree_prior(speciation_rate=0.4),
            "crown-conditioned-yule",
            {"speciation_rate": 0.4},
        ),
        (
            build_crown_conditioned_birth_death_tree_prior(
                speciation_rate=1.0,
                extinction_rate=0.25,
                sampling_fraction=0.8,
            ),
            "crown-conditioned-birth-death",
            {
                "speciation_rate": 1.0,
                "extinction_rate": 0.25,
                "sampling_fraction": 0.8,
            },
        ),
        (
            build_constant_population_coalescent_tree_prior(
                effective_population_size=1500.0
            ),
            "constant-population-coalescent",
            {"effective_population_size": 1500.0},
        ),
        (
            build_skyline_coalescent_tree_prior(
                [
                    SkylineCoalescentEpoch(
                        younger_boundary_age=0.0,
                        older_boundary_age=0.5,
                        effective_population_size=1000.0,
                    ),
                    SkylineCoalescentEpoch(
                        younger_boundary_age=0.5,
                        older_boundary_age=None,
                        effective_population_size=2500.0,
                    ),
                ]
            ),
            "skyline-coalescent",
            {
                "epoch_count": 2.0,
                "epoch_1_younger_boundary_age": 0.0,
                "epoch_1_older_boundary_age": 0.5,
                "epoch_1_effective_population_size": 1000.0,
                "epoch_2_younger_boundary_age": 0.5,
                "epoch_2_effective_population_size": 2500.0,
            },
        ),
    ],
)
def test_shared_prior_protocol_summarizes_supported_prior_models(
    prior_model: object,
    expected_family: str,
    expected_parameter_values: dict[str, float],
) -> None:
    assert supports_bayesian_prior_model_protocol(prior_model)

    summary = summarize_bayesian_prior_model(prior_model)

    assert summary == BayesianPriorModelSummary(
        family=expected_family,
        parameter_values=expected_parameter_values,
        parameter_count=len(expected_parameter_values),
    )


def test_shared_prior_protocol_supports_calibration_prior_definitions() -> None:
    calibration_prior = load_calibration_prior_definitions(
        fixture("trees", "calibration_prior_time_tree_6_taxa.nwk"),
        fixture("metadata", "calibration_priors_6_taxa.tsv"),
    )[0]

    assert supports_bayesian_prior_model_protocol(calibration_prior)
    assert summarize_bayesian_prior_model(calibration_prior) == BayesianPriorModelSummary(
        family="offset-exponential",
        parameter_values={
            "offset_age": 0.75,
            "exponential_mean": 1.0,
        },
        parameter_count=2,
    )


def test_shared_prior_protocol_rejects_objects_without_prior_fields() -> None:
    assert not supports_bayesian_prior_model_protocol(object())
