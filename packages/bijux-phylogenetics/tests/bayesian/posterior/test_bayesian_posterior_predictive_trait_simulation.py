from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.bayesian.brownian_continuous_trait import (
    build_brownian_continuous_trait_model_definition,
    build_brownian_continuous_trait_proposal_schedule,
    run_brownian_continuous_trait_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.continuous_trait_location_priors import (
    build_normal_continuous_trait_location_prior,
)
from bijux_phylogenetics.bayesian.continuous_trait_model_priors import (
    build_exponential_continuous_trait_scalar_prior,
)
from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    build_discrete_trait_mk_model_definition,
    build_discrete_trait_mk_proposal_schedule,
    run_discrete_trait_mk_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    build_exponential_discrete_trait_rate_prior,
)
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    build_ornstein_uhlenbeck_continuous_trait_model_definition,
    build_ornstein_uhlenbeck_continuous_trait_proposal_schedule,
    run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    build_posterior_predictive_simulation_definition,
    simulate_brownian_continuous_trait_posterior_predictive,
    simulate_discrete_trait_mk_posterior_predictive,
    simulate_ornstein_uhlenbeck_continuous_trait_posterior_predictive,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_discrete_trait_mk_posterior_predictive_emits_replicate_statistics() -> None:
    model_definition = build_discrete_trait_mk_model_definition(
        transition_model_name="symmetric",
        rate_prior=build_exponential_discrete_trait_rate_prior(rate=0.8),
    )
    run_report = run_discrete_trait_mk_metropolis_hastings(
        tree=_build_trait_tree(),
        tip_states=_build_tip_states(),
        model_definition=model_definition,
        proposal_schedule=build_discrete_trait_mk_proposal_schedule(
            model_definition=model_definition,
            rate_log_scale_standard_deviation=0.35,
        ),
        iteration_count=10,
        sample_every=1,
        seed=0,
    )

    report = simulate_discrete_trait_mk_posterior_predictive(
        run_report=run_report,
        tip_states=_build_tip_states(),
        definition=build_posterior_predictive_simulation_definition(
            replicate_count=3,
            seed=0,
        ),
    )

    assert report.model_name == "symmetric"
    assert len(report.observed_statistic_rows) == 2
    assert len(report.replicates) == 3
    assert len(report.replicate_statistic_rows) == 6
    assert len(report.statistic_summary_rows) == 2
    assert all(len(replicate.traits) == 4 for replicate in report.replicates)


def test_brownian_posterior_predictive_emits_replicate_statistics() -> None:
    model_definition = build_brownian_continuous_trait_model_definition(
        root_state_prior=build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=2.0,
        ),
        sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
    )
    run_report = run_brownian_continuous_trait_metropolis_hastings(
        tree=_load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
        tip_values={"A": 0.2, "B": 0.5, "C": 1.0, "D": 1.4},
        model_definition=model_definition,
        proposal_schedule=build_brownian_continuous_trait_proposal_schedule(
            model_definition=model_definition,
            root_state_move_weight=1.0,
            root_state_standard_deviation=0.25,
            sigma_squared_move_weight=1.0,
            sigma_squared_log_scale_standard_deviation=0.25,
        ),
        iteration_count=10,
        sample_every=1,
        seed=5,
    )

    report = simulate_brownian_continuous_trait_posterior_predictive(
        run_report=run_report,
        definition=build_posterior_predictive_simulation_definition(
            replicate_count=3,
            seed=0,
        ),
    )

    assert report.model_name == "brownian"
    assert len(report.observed_statistic_rows) == 3
    assert len(report.replicates) == 3
    assert len(report.replicate_statistic_rows) == 9
    assert len(report.statistic_summary_rows) == 3
    assert all(len(replicate.traits) == 4 for replicate in report.replicates)


def test_ou_posterior_predictive_emits_replicate_statistics() -> None:
    model_definition = build_ornstein_uhlenbeck_continuous_trait_model_definition(
        alpha_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
        optimum_prior=build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=2.0,
        ),
        sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
    )
    run_report = run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings(
        tree=_load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
        tip_values={"A": 0.2, "B": 0.5, "C": 1.0, "D": 1.4},
        model_definition=model_definition,
        proposal_schedule=build_ornstein_uhlenbeck_continuous_trait_proposal_schedule(
            model_definition=model_definition,
            alpha_move_weight=1.0,
            alpha_log_scale_standard_deviation=0.25,
            optimum_move_weight=1.0,
            optimum_standard_deviation=0.25,
            sigma_squared_move_weight=1.0,
            sigma_squared_log_scale_standard_deviation=0.25,
        ),
        iteration_count=10,
        sample_every=1,
        seed=5,
    )

    report = simulate_ornstein_uhlenbeck_continuous_trait_posterior_predictive(
        run_report=run_report,
        definition=build_posterior_predictive_simulation_definition(
            replicate_count=3,
            seed=0,
        ),
    )

    assert report.model_name == "ornstein-uhlenbeck"
    assert len(report.observed_statistic_rows) == 3
    assert len(report.replicates) == 3
    assert len(report.replicate_statistic_rows) == 9
    assert len(report.statistic_summary_rows) == 3
    assert all(len(replicate.traits) == 4 for replicate in report.replicates)


def _build_trait_tree() -> PhyloTree:
    tree = PhyloTree.from_newick("((A:0.2,B:0.2):0.3,(C:0.2,D:0.2):0.3);")
    tree.rooted = True
    return tree


def _build_tip_states() -> dict[str, str]:
    return {
        "A": "north",
        "B": "south",
        "C": "west",
        "D": "north",
    }


def _load_rooted_tree_fixture(name: str):
    tree = load_tree(fixture("trees", name))
    tree.rooted = True
    return tree
