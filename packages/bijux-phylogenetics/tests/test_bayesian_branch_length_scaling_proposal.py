from __future__ import annotations

import math
import random

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_branch_length_scaling_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def test_branch_length_scaling_proposal_keeps_all_branch_lengths_positive() -> None:
    current_state = _build_flat_branch_length_state(branch_lengths=(0.2, 0.3, 0.4, 0.5))
    rng = random.Random(17)

    for _ in range(100):
        proposal = propose_branch_length_scaling_move(
            current_state,
            rng,
            log_scale_standard_deviation=0.7,
        )
        assert proposal.is_valid is True
        assert proposal.invalid_reason is None
        assert len(proposal.changed_fields) == 1
        assert proposal.changed_fields[0].startswith("tree.branch_length:")
        assert proposal.proposed_tree is not None
        proposed_tree = proposal.proposed_tree
        assert all(
            child.branch_length is not None and child.branch_length > 0.0
            for _parent, child in proposed_tree.iter_edges()
        )


def test_branch_length_scaling_proposal_reports_hastings_correction() -> None:
    current_state = _build_flat_branch_length_state(branch_lengths=(0.2, 0.3, 0.4, 0.5))
    proposal = propose_branch_length_scaling_move(
        current_state,
        random.Random(9),
        log_scale_standard_deviation=0.5,
    )

    assert proposal.is_valid is True
    changed_branch_id = proposal.changed_fields[0].split("tree.branch_length:", 1)[1]
    current_branch_length = (
        current_state.tree.to_tree().node_by_id(changed_branch_id).branch_length
    )
    assert current_branch_length is not None
    proposed_tree = proposal.proposed_tree
    assert proposed_tree is not None
    proposed_branch_length = proposed_tree.node_by_id(changed_branch_id).branch_length
    assert proposed_branch_length is not None
    assert math.isclose(
        proposal.log_hastings_ratio,
        math.log(proposed_branch_length / current_branch_length),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_sampler_uses_branch_length_scaling_hastings_ratio_on_flat_target() -> None:
    initial_state = _build_flat_branch_length_state(branch_lengths=(0.2, 0.3, 0.4, 0.5))
    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: propose_branch_length_scaling_move(
            current_state,
            rng,
            log_scale_standard_deviation=0.5,
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=5,
        sample_every=1,
        seed=21,
    )

    for step_row in report.step_rows:
        assert step_row.proposal_valid is True
        assert step_row.log_acceptance_ratio is not None
        assert math.isclose(
            step_row.log_acceptance_ratio,
            step_row.log_hastings_ratio,
            rel_tol=0.0,
            abs_tol=1e-12,
        )


def test_branch_length_scaling_proposal_reports_invalid_tree_without_positive_edges() -> (
    None
):
    current_state = _build_flat_branch_length_state(branch_lengths=(0.0, 0.0, 0.0, 0.0))
    proposal = propose_branch_length_scaling_move(
        current_state,
        random.Random(5),
        log_scale_standard_deviation=0.5,
    )

    assert proposal.is_valid is False
    assert (
        proposal.invalid_reason == "tree has no strictly positive finite branch lengths"
    )
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def _build_flat_branch_length_state(
    *,
    branch_lengths: tuple[float, float, float, float],
) -> BayesianPhylogeneticState:
    left_tip_length, right_tip_length, left_internal_length, outgroup_length = (
        branch_lengths
    )
    return score_bayesian_phylogenetic_state(
        tree=PhyloTree(
            TreeNode(
                children=[
                    TreeNode(
                        children=[
                            TreeNode(name="A", branch_length=left_tip_length),
                            TreeNode(name="B", branch_length=right_tip_length),
                        ],
                        branch_length=left_internal_length,
                    ),
                    TreeNode(name="C", branch_length=outgroup_length),
                ]
            ),
            rooted=True,
        ),
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_zero_log_likelihood,
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
