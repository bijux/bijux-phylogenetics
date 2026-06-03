from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    build_discrete_trait_mk_model_definition,
    build_discrete_trait_mk_proposal_schedule,
    run_discrete_trait_mk_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    build_exponential_discrete_trait_rate_prior,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


@pytest.mark.parametrize(
    ("transition_model_name", "expected_parameter_count"),
    [
        ("equal-rates", 1),
        ("symmetric", 3),
        ("all-rates-different", 6),
    ],
)
def test_discrete_trait_mk_runner_samples_er_sym_and_ard_transition_rate_posteriors(
    transition_model_name: str,
    expected_parameter_count: int,
) -> None:
    model_definition = build_discrete_trait_mk_model_definition(
        transition_model_name=transition_model_name,
        rate_prior=build_exponential_discrete_trait_rate_prior(rate=0.8),
    )
    proposal_schedule = build_discrete_trait_mk_proposal_schedule(
        model_definition=model_definition,
        rate_log_scale_standard_deviation=0.35,
    )

    report = run_discrete_trait_mk_metropolis_hastings(
        tree=_build_trait_tree(),
        tip_states=_build_tip_states(),
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=12,
        sample_every=1,
        seed=0,
    )

    assert report.chain_report.accepted_count >= 1
    assert report.state_order == ["north", "south", "west"]
    assert len(report.posterior_rows) == len(report.chain_report.sampled_states)
    assert all(
        row.transition_model_name == transition_model_name
        for row in report.posterior_rows
    )
    assert all(
        len(row.rate_parameters) == expected_parameter_count
        for row in report.posterior_rows
    )
    assert all(
        len(row.prior_component_log_priors) == expected_parameter_count
        for row in report.posterior_rows
    )
    assert all(summary.node_state_summaries for summary in report.posterior_rows)
    assert all(
        any(
            node_summary.node_id == "A|B|C|D"
            for node_summary in row.node_state_summaries
        )
        for row in report.posterior_rows
    )
    assert all(
        state.model_parameters.categorical_parameters["discrete-trait-model"]
        == transition_model_name
        for state in report.chain_report.sampled_states
    )
    assert any(
        changed_field.startswith("vector_parameters.discrete-trait-rates.")
        for step_row in report.chain_report.step_rows
        for changed_field in step_row.proposal_changed_fields
    )


def test_discrete_trait_mk_runner_respects_fixed_root_state_in_node_state_summaries() -> (
    None
):
    model_definition = build_discrete_trait_mk_model_definition(
        transition_model_name="ER",
        rate_prior=build_exponential_discrete_trait_rate_prior(rate=0.8),
        root_prior_mode="fixed",
        fixed_root_state="north",
    )
    proposal_schedule = build_discrete_trait_mk_proposal_schedule(
        model_definition=model_definition,
        rate_log_scale_standard_deviation=0.3,
    )

    report = run_discrete_trait_mk_metropolis_hastings(
        tree=_build_trait_tree(),
        tip_states=_build_tip_states(),
        model_definition=model_definition,
        proposal_schedule=proposal_schedule,
        iteration_count=8,
        sample_every=1,
        seed=0,
    )

    root_summaries = [
        next(
            node_summary
            for node_summary in row.node_state_summaries
            if node_summary.node_id == "A|B|C|D"
        )
        for row in report.posterior_rows
    ]

    assert all(summary.most_likely_state == "north" for summary in root_summaries)
    assert all(
        summary.state_probabilities["north"] == 1.0 for summary in root_summaries
    )


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
