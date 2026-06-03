from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.adaptive_tuning import (
    AdaptiveMetropolisHastingsRunReport,
    build_adaptive_tuning_controller,
    run_adaptive_tuned_metropolis_hastings_sampler,
)
from bijux_phylogenetics.bayesian.metropolis_hastings import (
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


def test_adaptive_tuning_sampler_adjusts_scale_during_burnin_and_freezes_afterward() -> (
    None
):
    report = run_adaptive_tuned_metropolis_hastings_sampler(
        initial_state=_build_scored_normal_target_state(0.0),
        propose_state=_thresholded_scale_proposal,
        tuning_controller=build_adaptive_tuning_controller(
            proposal_name="thresholded-random-walk",
            scale_parameter_name="step-size",
            initial_scale=4.0,
            target_acceptance_rate=0.3,
            burnin_iteration_count=6,
            adaptation_window_size=2,
            decrease_factor=0.5,
            increase_factor=2.0,
            minimum_scale=0.25,
            maximum_scale=8.0,
        ),
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=8,
        sample_every=1,
        seed=7,
    )

    assert isinstance(report, AdaptiveMetropolisHastingsRunReport)
    assert report.freeze_iteration_index == 7
    assert report.burnin_iteration_count == 6
    assert report.burnin_sample_count == 7
    assert report.retained_sample_count == 2
    assert report.chain_report.accepted_count == 2
    assert report.chain_report.rejected_count == 6
    assert [row.action for row in report.tuning_report.window_rows] == [
        "decrease",
        "decrease",
        "increase",
        "frozen",
    ]
    assert [
        (row.scale_before_window, row.scale_after_window)
        for row in report.tuning_report.window_rows
    ] == [
        (4.0, 2.0),
        (2.0, 1.0),
        (1.0, 2.0),
        (2.0, 2.0),
    ]
    assert report.tuning_report.final_scale == 2.0
    assert all(
        row.within_burnin is False
        for row in report.tuning_report.window_rows
        if row.window_start_iteration >= report.freeze_iteration_index
    )


def test_adaptive_tuning_sampler_requires_retained_iterations_after_burnin() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="requires 'iteration_count' to exceed the tuning-controller burn-in",
    ):
        run_adaptive_tuned_metropolis_hastings_sampler(
            initial_state=_build_scored_normal_target_state(0.0),
            propose_state=_thresholded_scale_proposal,
            tuning_controller=build_adaptive_tuning_controller(
                proposal_name="thresholded-random-walk",
                scale_parameter_name="step-size",
                initial_scale=1.0,
                target_acceptance_rate=0.3,
                burnin_iteration_count=4,
                adaptation_window_size=2,
            ),
            update_prior_components=_normal_target_prior_components,
            update_log_likelihood=_zero_log_likelihood,
            iteration_count=4,
            sample_every=1,
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


def _thresholded_scale_proposal(
    current_state: BayesianPhylogeneticState,
    _rng,
    proposal_scale: float,
):
    if proposal_scale > 1.0:
        return build_metropolis_hastings_proposal(
            changed_fields=("scalar_parameters.x",),
            log_forward_density=10.0,
            log_reverse_density=0.0,
            is_valid=True,
            proposed_tree=current_state.tree.to_tree(),
            proposed_model_parameters=current_state.model_parameters,
        )
    return build_metropolis_hastings_proposal(
        changed_fields=("scalar_parameters.x",),
        log_forward_density=0.0,
        log_reverse_density=0.0,
        is_valid=True,
        proposed_tree=current_state.tree.to_tree(),
        proposed_model_parameters=current_state.model_parameters,
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
