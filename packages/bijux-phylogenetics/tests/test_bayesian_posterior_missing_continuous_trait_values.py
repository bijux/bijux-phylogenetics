from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    BrownianContinuousTraitRunReport,
    build_brownian_continuous_trait_model_definition,
    build_brownian_continuous_trait_proposal_schedule,
)
from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    build_normal_continuous_trait_location_prior,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_exponential_continuous_trait_scalar_prior,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import MetropolisHastingsRunReport
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    OrnsteinUhlenbeckContinuousTraitRunReport,
    build_ornstein_uhlenbeck_continuous_trait_model_definition,
    build_ornstein_uhlenbeck_continuous_trait_proposal_schedule,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    build_posterior_missing_continuous_trait_definition,
    summarize_brownian_continuous_trait_posterior_missing_values,
    summarize_continuous_trait_posterior_missing_values,
    summarize_ornstein_uhlenbeck_continuous_trait_posterior_missing_values,
)
from bijux_phylogenetics.bayesian.state import (
    build_bayesian_model_parameter_state,
    build_bayesian_phylogenetic_state,
    build_bayesian_prior_component_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_posterior_missing_continuous_trait_definition_requires_masked_taxon() -> None:
    with pytest.raises(PhylogeneticsError, match="at least one masked taxon"):
        build_posterior_missing_continuous_trait_definition(
            tip_values={"A": 1.0, "B": 2.0}
        )


def test_brownian_missing_value_summary_recovers_masked_known_value() -> None:
    report = summarize_brownian_continuous_trait_posterior_missing_values(
        run_report=_build_manual_brownian_run_report(),
        definition=build_posterior_missing_continuous_trait_definition(
            tip_values={"A": 1.2, "B": None, "C": 0.0}
        ),
    )

    assert report.sample_count == 1
    assert report.sampled_trait_models == ["brownian"]
    assert report.warnings == []

    summary_row = report.taxon_summary_rows[0]
    assert summary_row.taxon == "B"
    assert math.isclose(summary_row.posterior_mean, 1.06666666666667, rel_tol=1e-6)
    assert summary_row.posterior_hpd_95_lower < 1.1 < summary_row.posterior_hpd_95_upper
    assert summary_row.mean_conditional_standard_deviation > 0.0


def test_ou_missing_value_summary_recovers_masked_known_value() -> None:
    report = summarize_ornstein_uhlenbeck_continuous_trait_posterior_missing_values(
        run_report=_build_manual_ou_run_report(),
        definition=build_posterior_missing_continuous_trait_definition(
            tip_values={"A": 1.2, "B": float("nan"), "C": 0.0}
        ),
    )

    summary_row = report.taxon_summary_rows[0]
    assert summary_row.taxon == "B"
    assert 0.4 < summary_row.posterior_mean < 1.2
    assert summary_row.posterior_hpd_95_lower < 1.0 < summary_row.posterior_hpd_95_upper
    assert summary_row.mean_conditional_standard_deviation > 0.0


def test_generic_continuous_missing_value_summary_dispatches_ou_report() -> None:
    report = summarize_continuous_trait_posterior_missing_values(
        run_report=_build_manual_ou_run_report(),
        definition=build_posterior_missing_continuous_trait_definition(
            tip_values={"A": 1.2, "B": None, "C": 0.0}
        ),
    )

    assert report.sampled_trait_models == ["ornstein-uhlenbeck"]


def _build_manual_brownian_run_report() -> BrownianContinuousTraitRunReport:
    model_definition = build_brownian_continuous_trait_model_definition(
        root_state_prior=build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=2.0,
        ),
        sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
        initial_root_state=0.0,
        initial_sigma_squared=0.5,
    )
    proposal_schedule = build_brownian_continuous_trait_proposal_schedule(
        model_definition=model_definition,
        root_state_move_weight=1.0,
        root_state_standard_deviation=0.2,
        sigma_squared_move_weight=1.0,
        sigma_squared_log_scale_standard_deviation=0.2,
    )
    tree = PhyloTree.from_newick("((A:0.1,B:0.1):0.8,C:0.9);")
    tree.rooted = True
    sampled_state = build_bayesian_phylogenetic_state(
        tree=tree,
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"continuous-trait-model": "brownian"},
            scalar_parameters={"root-state": 0.0, "sigma-squared": 0.5},
        ),
        prior_components=[
            build_bayesian_prior_component_state(
                component_name="continuous-trait:root-state",
                family="normal",
                log_prior=-0.1,
            ),
            build_bayesian_prior_component_state(
                component_name="continuous-trait:sigma-squared",
                family="exponential",
                log_prior=-0.5,
            ),
        ],
        log_likelihood=-1.0,
    )
    return BrownianContinuousTraitRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        taxa=["A", "B", "C"],
        tip_values={"A": 1.2, "B": 1.1, "C": 0.0},
        chain_report=MetropolisHastingsRunReport(
            iteration_count=1,
            sample_every=1,
            seed=0,
            accepted_count=1,
            rejected_count=0,
            acceptance_rate=1.0,
            initial_state=sampled_state,
            final_state=sampled_state,
            sampled_states=[sampled_state],
            step_rows=[],
        ),
        posterior_rows=[],
        parameter_summaries=[],
    )


def _build_manual_ou_run_report() -> OrnsteinUhlenbeckContinuousTraitRunReport:
    model_definition = build_ornstein_uhlenbeck_continuous_trait_model_definition(
        alpha_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
        optimum_prior=build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=2.0,
        ),
        sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
        initial_alpha=0.8,
        initial_optimum=0.0,
        initial_sigma_squared=0.6,
    )
    proposal_schedule = build_ornstein_uhlenbeck_continuous_trait_proposal_schedule(
        model_definition=model_definition,
        alpha_move_weight=1.0,
        alpha_log_scale_standard_deviation=0.2,
        optimum_move_weight=1.0,
        optimum_standard_deviation=0.2,
        sigma_squared_move_weight=1.0,
        sigma_squared_log_scale_standard_deviation=0.2,
    )
    tree = PhyloTree.from_newick("((A:0.1,B:0.1):0.8,C:0.9);")
    tree.rooted = True
    sampled_state = build_bayesian_phylogenetic_state(
        tree=tree,
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"continuous-trait-model": "ornstein-uhlenbeck"},
            scalar_parameters={"alpha": 0.8, "optimum": 0.0, "sigma-squared": 0.6},
        ),
        prior_components=[
            build_bayesian_prior_component_state(
                component_name="continuous-trait:alpha",
                family="exponential",
                log_prior=-0.8,
            ),
            build_bayesian_prior_component_state(
                component_name="continuous-trait:optimum",
                family="normal",
                log_prior=-0.1,
            ),
            build_bayesian_prior_component_state(
                component_name="continuous-trait:sigma-squared",
                family="exponential",
                log_prior=-0.6,
            ),
        ],
        log_likelihood=-1.1,
    )
    return OrnsteinUhlenbeckContinuousTraitRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        taxa=["A", "B", "C"],
        tip_values={"A": 1.2, "B": 1.0, "C": 0.0},
        chain_report=MetropolisHastingsRunReport(
            iteration_count=1,
            sample_every=1,
            seed=0,
            accepted_count=1,
            rejected_count=0,
            acceptance_rate=1.0,
            initial_state=sampled_state,
            final_state=sampled_state,
            sampled_states=[sampled_state],
            step_rows=[],
        ),
        posterior_rows=[],
        parameter_summaries=[],
        identifiability_warnings=[],
    )
