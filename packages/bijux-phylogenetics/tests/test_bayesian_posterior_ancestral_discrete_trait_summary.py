from __future__ import annotations

from dataclasses import replace
import math

import pytest

from bijux_phylogenetics.bayesian.discrete_trait_mk import (
    DiscreteTraitMkNodeStateSummary,
    build_discrete_trait_mk_model_definition,
    build_discrete_trait_mk_proposal_schedule,
    run_discrete_trait_mk_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.discrete_trait_rate_priors import (
    build_exponential_discrete_trait_rate_prior,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_traits import (
    summarize_discrete_trait_mk_posterior_ancestral_states,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_discrete_trait_posterior_summary_emits_state_probabilities_and_entropy() -> (
    None
):
    model_definition = build_discrete_trait_mk_model_definition(
        transition_model_name="SYM",
        rate_prior=build_exponential_discrete_trait_rate_prior(rate=0.8),
    )
    report = run_discrete_trait_mk_metropolis_hastings(
        tree=_build_trait_tree(),
        tip_states=_build_tip_states(),
        model_definition=model_definition,
        proposal_schedule=build_discrete_trait_mk_proposal_schedule(
            model_definition=model_definition,
            rate_log_scale_standard_deviation=0.35,
        ),
        iteration_count=12,
        sample_every=1,
        seed=3,
    )

    summary = summarize_discrete_trait_mk_posterior_ancestral_states(report)

    assert summary.sample_count == len(report.posterior_rows)
    assert summary.distinct_topology_count == 1
    assert summary.sampled_transition_models == ["symmetric"]
    assert summary.state_order == ["north", "south", "west"]
    assert summary.tree_uncertainty_policy == "fixed-topology-posterior-aggregation"
    assert summary.warnings == []

    root_state_rows = [
        row
        for row in summary.state_probability_rows
        if row.descendant_taxa == ["A", "B", "C", "D"]
    ]
    assert [row.state for row in root_state_rows] == ["north", "south", "west"]
    assert math.isclose(
        sum(row.conditional_posterior_probability for row in root_state_rows),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-9,
    )
    assert math.isclose(
        sum(row.marginal_posterior_probability for row in root_state_rows),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-9,
    )
    assert all(row.clade_posterior_probability == 1.0 for row in root_state_rows)
    assert all(
        row.supporting_sample_count == summary.sample_count for row in root_state_rows
    )

    root_summary = next(
        row
        for row in summary.node_summary_rows
        if row.descendant_taxa == ["A", "B", "C", "D"]
    )
    assert root_summary.most_likely_state in {"north", "south", "west"}
    assert root_summary.clade_posterior_probability == 1.0
    assert 0.0 <= root_summary.max_conditional_posterior_probability <= 1.0
    assert 0.0 <= root_summary.max_marginal_posterior_probability <= 1.0
    assert 0.0 <= root_summary.conditional_posterior_entropy <= math.log(3.0)


def test_discrete_trait_posterior_summary_rejects_missing_state_probability() -> None:
    model_definition = build_discrete_trait_mk_model_definition(
        transition_model_name="ER",
        rate_prior=build_exponential_discrete_trait_rate_prior(rate=0.8),
    )
    report = run_discrete_trait_mk_metropolis_hastings(
        tree=_build_trait_tree(),
        tip_states=_build_tip_states(),
        model_definition=model_definition,
        proposal_schedule=build_discrete_trait_mk_proposal_schedule(
            model_definition=model_definition,
            rate_log_scale_standard_deviation=0.3,
        ),
        iteration_count=4,
        sample_every=1,
        seed=1,
    )
    first_row = report.posterior_rows[0]
    broken_node_summary = DiscreteTraitMkNodeStateSummary(
        node_id=first_row.node_state_summaries[0].node_id,
        node_name=first_row.node_state_summaries[0].node_name,
        descendant_taxa=list(first_row.node_state_summaries[0].descendant_taxa),
        most_likely_state=first_row.node_state_summaries[0].most_likely_state,
        state_probabilities={"north": 0.7, "south": 0.3},
    )
    broken_row = replace(
        first_row,
        node_state_summaries=[broken_node_summary, *first_row.node_state_summaries[1:]],
    )
    broken_report = replace(
        report,
        posterior_rows=[broken_row, *report.posterior_rows[1:]],
    )

    with pytest.raises(PhylogeneticsError, match="report every modeled discrete state"):
        summarize_discrete_trait_mk_posterior_ancestral_states(broken_report)


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
