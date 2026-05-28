from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    MetropolisHastingsProposal,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_metropolis_hastings_recovers_standard_normal_target() -> None:
    initial_state = _build_scored_normal_target_state(0.0)

    def propose_state(
        current_state: BayesianPhylogeneticState,
        rng,
    ) -> MetropolisHastingsProposal:
        current_x = current_state.model_parameters.scalar_parameters["x"]
        proposed_x = current_x + rng.gauss(0.0, 1.0)
        return MetropolisHastingsProposal(
            proposed_tree=current_state.tree.to_tree(),
            proposed_model_parameters=build_bayesian_model_parameter_state(
                scalar_parameters={"x": proposed_x}
            ),
            log_hastings_ratio=0.0,
        )

    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=propose_state,
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=12000,
        sample_every=5,
        seed=23,
    )

    sampled_x = [
        state.model_parameters.scalar_parameters["x"]
        for state in report.sampled_states[401:]
    ]
    sample_mean = sum(sampled_x) / len(sampled_x)
    sample_variance = sum((value - sample_mean) ** 2 for value in sampled_x) / len(
        sampled_x
    )

    assert abs(sample_mean) < 0.12
    assert abs(sample_variance - 1.0) < 0.18
    assert report.accepted_count > 0
    assert report.rejected_count > 0
    assert 0.2 < report.acceptance_rate < 0.8


def test_metropolis_hastings_uses_hastings_ratio_in_acceptance_decision() -> None:
    initial_state = _build_scored_normal_target_state(0.0)

    def zero_shift_proposal(
        current_state: BayesianPhylogeneticState,
        _rng,
    ) -> MetropolisHastingsProposal:
        return MetropolisHastingsProposal(
            proposed_tree=current_state.tree.to_tree(),
            proposed_model_parameters=current_state.model_parameters,
            log_hastings_ratio=0.0,
        )

    def penalized_zero_shift_proposal(
        current_state: BayesianPhylogeneticState,
        _rng,
    ) -> MetropolisHastingsProposal:
        return MetropolisHastingsProposal(
            proposed_tree=current_state.tree.to_tree(),
            proposed_model_parameters=current_state.model_parameters,
            log_hastings_ratio=-10.0,
        )

    accepted_report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=zero_shift_proposal,
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=5,
        sample_every=1,
        seed=7,
    )
    rejected_report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=penalized_zero_shift_proposal,
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=5,
        sample_every=1,
        seed=7,
    )

    assert all(step.accepted for step in accepted_report.step_rows)
    assert all(not step.accepted for step in rejected_report.step_rows)
    assert all(
        math.isclose(step.log_hastings_ratio, -10.0, rel_tol=0.0, abs_tol=1e-12)
        for step in rejected_report.step_rows
    )


def test_metropolis_hastings_requires_positive_iteration_count() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="requires 'iteration_count' to be positive",
    ):
        run_metropolis_hastings_sampler(
            initial_state=_build_scored_normal_target_state(0.0),
            propose_state=lambda current_state, rng: MetropolisHastingsProposal(
                proposed_tree=current_state.tree.to_tree(),
                proposed_model_parameters=current_state.model_parameters,
                log_hastings_ratio=0.0,
            ),
            update_prior_components=_normal_target_prior_components,
            update_log_likelihood=_zero_log_likelihood,
            iteration_count=0,
            seed=0,
        )


def test_metropolis_hastings_requires_finite_hastings_ratio() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="requires 'log_hastings_ratio' to be finite",
    ):
        run_metropolis_hastings_sampler(
            initial_state=_build_scored_normal_target_state(0.0),
            propose_state=lambda current_state, rng: MetropolisHastingsProposal(
                proposed_tree=current_state.tree.to_tree(),
                proposed_model_parameters=current_state.model_parameters,
                log_hastings_ratio=math.inf,
            ),
            update_prior_components=_normal_target_prior_components,
            update_log_likelihood=_zero_log_likelihood,
            iteration_count=1,
            seed=0,
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
