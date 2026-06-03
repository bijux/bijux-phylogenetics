from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    build_metropolis_hastings_proposal,
    deserialize_metropolis_hastings_checkpoint,
    deserialize_metropolis_hastings_checkpoint_json,
    run_checkpointed_metropolis_hastings_sampler,
    serialize_metropolis_hastings_checkpoint,
    serialize_metropolis_hastings_checkpoint_json,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_metropolis_hastings_checkpoint_round_trip_preserves_trace_state() -> None:
    checkpointed_report = run_checkpointed_metropolis_hastings_sampler(
        initial_state=_build_scored_normal_target_state(0.0),
        propose_state=_propose_normal_random_walk_state,
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=9,
        sample_every=2,
        seed=19,
        checkpoint_iteration_indexes=(3, 6, 9),
        stop_after_iteration_index=6,
    )

    assert checkpointed_report.completed is False
    assert checkpointed_report.completed_iteration_count == 6
    assert len(checkpointed_report.checkpoints) == 2

    checkpoint = checkpointed_report.checkpoints[-1]
    payload = serialize_metropolis_hastings_checkpoint(checkpoint)
    round_tripped_from_mapping = deserialize_metropolis_hastings_checkpoint(payload)
    round_tripped_from_json = deserialize_metropolis_hastings_checkpoint_json(
        serialize_metropolis_hastings_checkpoint_json(checkpoint)
    )

    assert round_tripped_from_mapping == checkpoint
    assert round_tripped_from_json == checkpoint
    assert checkpoint.current_state == checkpointed_report.current_state
    assert checkpoint.completed_iteration_count == 6
    assert len(checkpoint.sampled_states) == 4


def test_metropolis_hastings_checkpoint_rejects_inconsistent_trace_counts() -> None:
    checkpointed_report = run_checkpointed_metropolis_hastings_sampler(
        initial_state=_build_scored_normal_target_state(0.0),
        propose_state=_propose_normal_random_walk_state,
        update_prior_components=_normal_target_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=7,
        sample_every=2,
        seed=5,
        checkpoint_iteration_indexes=(4,),
        stop_after_iteration_index=4,
    )

    broken_payload = serialize_metropolis_hastings_checkpoint(
        checkpointed_report.checkpoints[-1]
    )
    broken_payload["accepted_count"] = 99

    with pytest.raises(
        PhylogeneticsError,
        match="accepted_count",
    ):
        deserialize_metropolis_hastings_checkpoint(broken_payload)


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
    from bijux_phylogenetics.bayesian.metropolis_hastings import (
        score_bayesian_phylogenetic_state,
    )

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
