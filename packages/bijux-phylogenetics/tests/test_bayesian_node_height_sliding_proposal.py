from __future__ import annotations

import math
import random

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_node_height_sliding_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def test_node_height_sliding_proposal_preserves_tip_depths() -> None:
    current_state = _build_ultrametric_node_height_state()
    proposal = propose_node_height_sliding_move(
        current_state,
        random.Random(11),
        height_slide_standard_deviation=0.15,
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert len(proposal.changed_fields) == 1
    assert proposal.changed_fields[0].startswith("tree.node_height:")
    proposed_tree = proposal.proposed_tree
    assert proposed_tree is not None
    current_tree = current_state.tree.to_tree()
    assert proposed_tree.root_to_tip_lengths() == current_tree.root_to_tip_lengths()
    current_branch_lengths = [
        child.branch_length for _parent, child in current_tree.iter_edges()
    ]
    proposed_branch_lengths = [
        child.branch_length for _parent, child in proposed_tree.iter_edges()
    ]
    assert proposed_branch_lengths[3] == current_branch_lengths[3]
    assert all(
        branch_length is not None and branch_length > 0.0
        for branch_length in proposed_branch_lengths
    )
    assert proposed_branch_lengths[:3] != current_branch_lengths[:3]


def test_node_height_sliding_proposal_is_symmetric_on_simple_ultrametric_tree() -> None:
    current_state = _build_ultrametric_node_height_state()
    proposal = propose_node_height_sliding_move(
        current_state,
        random.Random(19),
        height_slide_standard_deviation=0.15,
    )

    assert proposal.is_valid is True
    assert math.isclose(
        proposal.log_hastings_ratio,
        0.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_invalid_node_height_order_proposal_is_rejected_before_scoring() -> None:
    scored_callback_count = 0

    def update_prior_components(
        state: BayesianPhylogeneticState,
    ) -> list[BayesianPriorComponentState]:
        nonlocal scored_callback_count
        scored_callback_count += 1
        return _zero_prior_components(state)

    def update_log_likelihood(state: BayesianPhylogeneticState) -> float:
        nonlocal scored_callback_count
        scored_callback_count += 1
        return _zero_log_likelihood(state)

    initial_state = score_bayesian_phylogenetic_state(
        tree=_build_ultrametric_node_height_tree(),
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=update_prior_components,
        update_log_likelihood=update_log_likelihood,
    )
    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: propose_node_height_sliding_move(
            current_state,
            rng,
            height_slide_standard_deviation=0.5,
        ),
        update_prior_components=update_prior_components,
        update_log_likelihood=update_log_likelihood,
        iteration_count=1,
        sample_every=1,
        seed=1,
    )

    assert scored_callback_count == 2
    assert report.accepted_count == 0
    assert report.rejected_count == 1
    assert report.step_rows[0].proposal_valid is False
    assert report.step_rows[0].proposal_invalid_reason == (
        "node-height sliding would violate ancestor-descendant age order"
    )
    assert report.step_rows[0].proposed_posterior_log_score is None
    assert report.step_rows[0].log_acceptance_ratio is None


def test_node_height_sliding_proposal_rejects_non_ultrametric_tree() -> None:
    current_state = score_bayesian_phylogenetic_state(
        tree=PhyloTree(
            TreeNode(
                children=[
                    TreeNode(
                        children=[
                            TreeNode(name="A", branch_length=0.4),
                            TreeNode(name="B", branch_length=0.5),
                        ],
                        branch_length=0.6,
                    ),
                    TreeNode(name="C", branch_length=1.0),
                ]
            ),
            rooted=True,
        ),
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )
    proposal = propose_node_height_sliding_move(
        current_state,
        random.Random(5),
        height_slide_standard_deviation=0.15,
    )

    assert proposal.is_valid is False
    assert (
        proposal.invalid_reason
        == "node-height sliding requires one rooted ultrametric tree"
    )
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def _build_ultrametric_node_height_state() -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=_build_ultrametric_node_height_tree(),
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )


def _build_ultrametric_node_height_tree() -> PhyloTree:
    return PhyloTree(
        TreeNode(
            children=[
                TreeNode(
                    children=[
                        TreeNode(name="A", branch_length=0.4),
                        TreeNode(name="B", branch_length=0.4),
                    ],
                    branch_length=0.6,
                ),
                TreeNode(name="C", branch_length=1.0),
            ]
        ),
        rooted=True,
    )


def _zero_prior_components(
    _state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    return [
        BayesianPriorComponentState(
            component_name="flat-prior",
            family="constant",
            log_prior=0.0,
        )
    ]


def _zero_log_likelihood(_state: BayesianPhylogeneticState) -> float:
    return 0.0
