from __future__ import annotations

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    build_metropolis_hastings_proposal,
    resume_metropolis_hastings_sampler,
    run_checkpointed_metropolis_hastings_sampler,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def test_resumed_metropolis_hastings_chain_matches_uninterrupted_chain() -> None:
    checkpoint_schedule = (4, 8, 12)
    full_run = run_metropolis_hastings_sampler(
        initial_state=_build_scored_normal_target_state(0.0),
        propose_state=_propose_normal_random_walk_state,
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=12,
        sample_every=3,
        seed=23,
    )
    interrupted_run = run_checkpointed_metropolis_hastings_sampler(
        initial_state=_build_scored_normal_target_state(0.0),
        propose_state=_propose_normal_random_walk_state,
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=12,
        sample_every=3,
        seed=23,
        checkpoint_iteration_indexes=checkpoint_schedule,
        stop_after_iteration_index=8,
    )

    resumed_run = resume_metropolis_hastings_sampler(
        checkpoint=interrupted_run.checkpoints[-1],
        propose_state=_propose_normal_random_walk_state,
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        checkpoint_iteration_indexes=checkpoint_schedule,
    )
    resumed_chain = resumed_run.to_chain_report()

    assert interrupted_run.resumed is False
    assert interrupted_run.completed is False
    assert interrupted_run.completed_iteration_count == 8
    assert resumed_run.resumed is True
    assert resumed_run.completed is True
    assert [
        checkpoint.completed_iteration_count for checkpoint in resumed_run.checkpoints
    ] == [12]
    assert resumed_chain == full_run


def _propose_normal_random_walk_state(
    current_state: BayesianPhylogeneticState,
    rng,
):
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
