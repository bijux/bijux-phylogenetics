from __future__ import annotations

import math
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
from bijux_phylogenetics.bayesian.ornstein_uhlenbeck_continuous_trait import (
    build_ornstein_uhlenbeck_continuous_trait_model_definition,
    build_ornstein_uhlenbeck_continuous_trait_proposal_schedule,
    run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_traits import (
    summarize_brownian_continuous_trait_posterior_ancestral_states,
    summarize_ornstein_uhlenbeck_continuous_trait_posterior_ancestral_states,
)
from bijux_phylogenetics.io.trees import load_tree

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_brownian_posterior_ancestral_summary_emits_root_and_internal_hpd_rows() -> (
    None
):
    model_definition = build_brownian_continuous_trait_model_definition(
        root_state_prior=build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=2.0,
        ),
        sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
    )
    report = run_brownian_continuous_trait_metropolis_hastings(
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
        iteration_count=20,
        sample_every=1,
        seed=5,
    )

    summary = summarize_brownian_continuous_trait_posterior_ancestral_states(report)

    assert summary.sample_count == len(report.posterior_rows)
    assert summary.distinct_topology_count == 1
    assert summary.sampled_trait_models == ["brownian"]
    assert summary.tree_uncertainty_policy == "fixed-topology-posterior-aggregation"
    assert summary.warnings == []

    root_row = next(
        row
        for row in summary.node_summary_rows
        if row.descendant_taxa == ["A", "B", "C", "D"]
    )
    internal_row = next(
        row for row in summary.node_summary_rows if row.descendant_taxa == ["A", "B"]
    )
    root_parameter_summary = next(
        parameter_summary
        for parameter_summary in report.parameter_summaries
        if parameter_summary.parameter_name == "root-state"
    )

    assert root_row.clade_posterior_probability == 1.0
    assert root_row.mean_conditional_standard_deviation == 0.0
    assert math.isclose(
        root_row.conditional_posterior_mean,
        root_parameter_summary.posterior_mean,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert root_row.conditional_hpd_95_lower <= root_row.conditional_hpd_95_upper
    assert internal_row.clade_posterior_probability == 1.0
    assert internal_row.mean_conditional_standard_deviation > 0.0
    assert (
        internal_row.conditional_hpd_95_lower <= internal_row.conditional_hpd_95_upper
    )


def test_ou_posterior_ancestral_summary_emits_optimum_root_and_internal_hpd_rows() -> (
    None
):
    model_definition = build_ornstein_uhlenbeck_continuous_trait_model_definition(
        alpha_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
        optimum_prior=build_normal_continuous_trait_location_prior(
            mean=0.0,
            standard_deviation=2.0,
        ),
        sigma_squared_prior=build_exponential_continuous_trait_scalar_prior(rate=1.0),
    )
    report = run_ornstein_uhlenbeck_continuous_trait_metropolis_hastings(
        tree=_load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk"),
        tip_values={"A": 0.2, "B": 0.5, "C": 1.0, "D": 1.4},
        model_definition=model_definition,
        proposal_schedule=build_ornstein_uhlenbeck_continuous_trait_proposal_schedule(
            model_definition=model_definition,
            alpha_move_weight=1.0,
            alpha_log_scale_standard_deviation=0.4,
            optimum_move_weight=1.0,
            optimum_standard_deviation=0.25,
            sigma_squared_move_weight=1.0,
            sigma_squared_log_scale_standard_deviation=0.4,
        ),
        iteration_count=24,
        sample_every=1,
        seed=7,
    )

    summary = summarize_ornstein_uhlenbeck_continuous_trait_posterior_ancestral_states(
        report
    )

    assert summary.sample_count == len(report.posterior_rows)
    assert summary.distinct_topology_count == 1
    assert summary.sampled_trait_models == ["ornstein-uhlenbeck"]
    assert summary.tree_uncertainty_policy == "fixed-topology-posterior-aggregation"
    assert summary.warnings == []

    root_row = next(
        row
        for row in summary.node_summary_rows
        if row.descendant_taxa == ["A", "B", "C", "D"]
    )
    internal_row = next(
        row for row in summary.node_summary_rows if row.descendant_taxa == ["A", "B"]
    )
    optimum_parameter_summary = next(
        parameter_summary
        for parameter_summary in report.parameter_summaries
        if parameter_summary.parameter_name == "optimum"
    )

    assert root_row.clade_posterior_probability == 1.0
    assert root_row.mean_conditional_standard_deviation == 0.0
    assert math.isclose(
        root_row.conditional_posterior_mean,
        optimum_parameter_summary.posterior_mean,
        rel_tol=0.0,
        abs_tol=1e-6,
    )
    assert root_row.conditional_hpd_95_lower <= root_row.conditional_hpd_95_upper
    assert internal_row.clade_posterior_probability == 1.0
    assert internal_row.mean_conditional_standard_deviation > 0.0
    assert (
        internal_row.conditional_hpd_95_lower <= internal_row.conditional_hpd_95_upper
    )


def _load_rooted_tree_fixture(name: str):
    tree = load_tree(fixture("trees", name))
    tree.rooted = True
    return tree
