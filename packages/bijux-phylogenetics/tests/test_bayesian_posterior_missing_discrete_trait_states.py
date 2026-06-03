from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    DiscreteTraitMkRunReport,
    build_discrete_trait_mk_model_definition,
    build_discrete_trait_mk_proposal_schedule,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    build_exponential_discrete_trait_rate_prior,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import MetropolisHastingsRunReport
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    build_posterior_missing_discrete_trait_definition,
    summarize_discrete_trait_mk_posterior_missing_states,
)
from bijux_phylogenetics.bayesian.state import (
    build_bayesian_model_parameter_state,
    build_bayesian_phylogenetic_state,
    build_bayesian_prior_component_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_posterior_missing_discrete_trait_definition_requires_masked_taxon() -> None:
    with pytest.raises(PhylogeneticsError, match="at least one masked taxon"):
        build_posterior_missing_discrete_trait_definition(
            tip_states={"A": "present", "B": "absent"}
        )


def test_discrete_trait_missing_state_summary_recovers_masked_known_state() -> None:
    report = summarize_discrete_trait_mk_posterior_missing_states(
        run_report=_build_manual_discrete_trait_run_report(),
        definition=build_posterior_missing_discrete_trait_definition(
            tip_states={"A": "present", "B": "?", "C": "absent"},
        ),
    )

    assert report.sample_count == 1
    assert report.sampled_transition_models == ["equal-rates"]
    assert report.warnings == []

    summary_row = report.taxon_summary_rows[0]
    assert summary_row.taxon == "B"
    assert summary_row.observed_symbol == "?"
    assert summary_row.most_likely_state == "present"
    assert 0.5 < summary_row.max_posterior_probability < 1.0
    assert summary_row.posterior_entropy > 0.0

    probability_by_state = {
        row.state: row.posterior_probability for row in report.state_probability_rows
    }
    assert probability_by_state["present"] > probability_by_state["absent"]


def _build_manual_discrete_trait_run_report() -> DiscreteTraitMkRunReport:
    model_definition = build_discrete_trait_mk_model_definition(
        transition_model_name="ER",
        rate_prior=build_exponential_discrete_trait_rate_prior(rate=1.0),
        root_prior_mode="equal",
        initial_rate=0.8,
    )
    proposal_schedule = build_discrete_trait_mk_proposal_schedule(
        model_definition=model_definition,
        rate_log_scale_standard_deviation=0.3,
    )
    tree = PhyloTree.from_newick("((A:0.05,B:0.05):0.8,C:0.85);")
    tree.rooted = True
    sampled_state = build_bayesian_phylogenetic_state(
        tree=tree,
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={
                "discrete-trait-model": "equal-rates",
                "root-prior-mode": "equal",
            },
            vector_parameters={"discrete-trait-rates": {"shared-rate": 0.8}},
        ),
        prior_components=[
            build_bayesian_prior_component_state(
                component_name="discrete-trait-rate:shared-rate",
                family="exponential",
                log_prior=-0.8,
            )
        ],
        log_likelihood=-1.2,
    )
    return DiscreteTraitMkRunReport(
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        state_order=["absent", "present"],
        taxa=["A", "B", "C"],
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
    )
