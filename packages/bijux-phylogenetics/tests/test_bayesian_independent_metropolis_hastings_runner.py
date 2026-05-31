from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.independent_chains import (
    build_independent_metropolis_hastings_chain_definition,
    run_independent_metropolis_hastings_chains,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsProposal,
    build_metropolis_hastings_proposal,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_independent_metropolis_hastings_runner_records_per_chain_outputs() -> None:
    report = run_independent_metropolis_hastings_chains(
        initial_state=_build_scored_normal_target_state(0.0),
        propose_state=_propose_normal_random_walk_state,
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        chain_definitions=(
            build_independent_metropolis_hastings_chain_definition(
                chain_name="chain_alpha",
                seed=11,
            ),
            build_independent_metropolis_hastings_chain_definition(
                chain_name="chain_beta",
                seed=19,
            ),
        ),
        iteration_count=8,
        sample_every=2,
    )

    assert [chain_report.chain_name for chain_report in report.chain_reports] == [
        "chain_alpha",
        "chain_beta",
    ]
    assert [chain_report.seed for chain_report in report.chain_reports] == [11, 19]
    assert all(
        len(chain_report.step_rows) == 8 for chain_report in report.chain_reports
    )
    assert all(
        len(chain_report.sampled_states) == 5 for chain_report in report.chain_reports
    )
    assert all(
        chain_report.acceptance_rate == chain_report.chain_report.acceptance_rate
        for chain_report in report.chain_reports
    )
    assert (
        report.chain_reports[0].trace_fingerprint
        != report.chain_reports[1].trace_fingerprint
    )
    assert report.diagnostics.chain_count == 2
    assert len(report.diagnostics.comparison_rows) == 1


def test_independent_metropolis_hastings_runner_requires_distinct_seeds() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="distinct seeds",
    ):
        run_independent_metropolis_hastings_chains(
            initial_state=_build_scored_normal_target_state(0.0),
            propose_state=_propose_normal_random_walk_state,
            update_prior_components=_normal_target_prior_components,
            update_log_likelihood=_zero_log_likelihood,
            chain_definitions=(
                build_independent_metropolis_hastings_chain_definition(
                    chain_name="chain_alpha",
                    seed=11,
                ),
                build_independent_metropolis_hastings_chain_definition(
                    chain_name="chain_beta",
                    seed=11,
                ),
            ),
            iteration_count=4,
        )


def _propose_normal_random_walk_state(
    current_state: BayesianPhylogeneticState,
    rng,
) -> MetropolisHastingsProposal:
    current_x = current_state.model_parameters.scalar_parameters["x"]
    proposed_x = current_x + rng.gauss(0.0, 1.0)
    return build_metropolis_hastings_proposal(
        changed_fields=("scalar_parameters.x",),
        log_forward_density=0.0,
        log_reverse_density=0.0,
        is_valid=True,
        proposed_tree=current_state.tree.to_tree(),
        proposed_model_parameters=build_bayesian_model_parameter_state(
            scalar_parameters={"x": proposed_x}
        ),
    )


def _build_scored_normal_target_state(x_value: float) -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=PhyloTree(
            TreeNode(
                children=[
                    TreeNode(name="A", branch_length=0.1),
                    TreeNode(name="B", branch_length=0.2),
                ]
            ),
            rooted=True,
        ),
        model_parameters=build_bayesian_model_parameter_state(
            scalar_parameters={"x": x_value}
        ),
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )


def _normal_target_prior_components(
    state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    x_value = state.model_parameters.scalar_parameters["x"]
    return [
        BayesianPriorComponentState(
            component_name="normal-target",
            family="gaussian",
            log_prior=-0.5 * (x_value**2),
            parameter_values={"mean": 0.0, "variance": 1.0},
        )
    ]


def _zero_log_likelihood(_state: BayesianPhylogeneticState) -> float:
    return 0.0
