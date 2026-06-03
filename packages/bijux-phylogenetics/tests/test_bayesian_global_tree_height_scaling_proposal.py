from __future__ import annotations

import math
import random

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_global_tree_height_scaling_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def test_global_tree_height_scaling_proposal_scales_all_branches_coherently() -> None:
    current_state = _build_positive_branch_state(branch_lengths=(0.2, 0.3, 0.4, 0.5))
    proposal = propose_global_tree_height_scaling_move(
        current_state,
        random.Random(11),
        log_scale_standard_deviation=0.4,
    )

    assert proposal.is_valid is True
    assert proposal.changed_fields == ("tree.branch_lengths",)
    proposed_tree = proposal.proposed_tree
    assert proposed_tree is not None
    current_branch_lengths = [
        float(child.branch_length)
        for _parent, child in current_state.tree.to_tree().iter_edges()
        if child.branch_length is not None
    ]
    proposed_branch_lengths = [
        float(child.branch_length)
        for _parent, child in proposed_tree.iter_edges()
        if child.branch_length is not None
    ]
    scale_factors = [
        proposed / current
        for current, proposed in zip(
            current_branch_lengths, proposed_branch_lengths, strict=True
        )
    ]
    first_scale_factor = scale_factors[0]
    assert all(branch_length > 0.0 for branch_length in proposed_branch_lengths)
    assert all(
        math.isclose(scale_factor, first_scale_factor, rel_tol=0.0, abs_tol=1e-12)
        for scale_factor in scale_factors
    )


def test_global_tree_height_scaling_proposal_reports_hastings_correction() -> None:
    current_state = _build_positive_branch_state(branch_lengths=(0.2, 0.3, 0.4, 0.5))
    proposal = propose_global_tree_height_scaling_move(
        current_state,
        random.Random(19),
        log_scale_standard_deviation=0.5,
    )

    assert proposal.is_valid is True
    proposed_tree = proposal.proposed_tree
    assert proposed_tree is not None
    current_branch_lengths = [
        float(child.branch_length)
        for _parent, child in current_state.tree.to_tree().iter_edges()
        if child.branch_length is not None
    ]
    proposed_branch_lengths = [
        float(child.branch_length)
        for _parent, child in proposed_tree.iter_edges()
        if child.branch_length is not None
    ]
    scale_factor = proposed_branch_lengths[0] / current_branch_lengths[0]
    assert math.isclose(
        proposal.log_hastings_ratio,
        len(current_branch_lengths) * math.log(scale_factor),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_sampler_uses_global_tree_height_scaling_hastings_ratio_on_flat_target() -> (
    None
):
    initial_state = _build_positive_branch_state(branch_lengths=(0.2, 0.3, 0.4, 0.5))
    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: (
            propose_global_tree_height_scaling_move(
                current_state,
                rng,
                log_scale_standard_deviation=0.5,
            )
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=5,
        sample_every=1,
        seed=29,
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


def test_global_tree_height_scaling_proposal_rejects_non_positive_branch_lengths() -> (
    None
):
    current_state = _build_positive_branch_state(branch_lengths=(0.0, 0.3, 0.4, 0.5))
    proposal = propose_global_tree_height_scaling_move(
        current_state,
        random.Random(7),
        log_scale_standard_deviation=0.5,
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "global tree-height scaling requires strictly positive finite branch lengths on every edge"
    )
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def _build_positive_branch_state(
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
