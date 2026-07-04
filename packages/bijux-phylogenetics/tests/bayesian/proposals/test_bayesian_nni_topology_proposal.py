from __future__ import annotations

from random import Random

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_nni_topology_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.phylo.topology.clades import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def test_nni_topology_proposal_changes_rooted_topology_without_changing_taxa() -> None:
    current_state = _build_scored_nni_state(
        update_prior_components=_zero_prior_components
    )
    current_tree = current_state.tree.to_tree()
    current_topology_fingerprint = rooted_topology_fingerprint(current_tree)

    proposal = propose_nni_topology_move(current_state, Random(1))

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.changed_fields == (
        "tree.topology:nni:root:clade:A|B|C:root:clade:A|B|C/clade:A|B:root:clade:A|B|C/clade:A|B/taxon:A",
    )
    proposed_tree = proposal.proposed_tree
    assert proposed_tree is not None
    assert set(proposed_tree.tip_names) == set(current_tree.tip_names)
    assert all(
        len(node.children) == 2
        for node in proposed_tree.iter_internal_nodes(order="preorder")
    )
    assert rooted_topology_fingerprint(proposed_tree) != current_topology_fingerprint


def test_nni_topology_sampler_scores_and_accepts_rejects_moves() -> None:
    initial_state = _build_scored_nni_state(
        update_prior_components=_topology_preference_prior_components
    )

    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=propose_nni_topology_move,
        update_prior_components=_topology_preference_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=6,
        sample_every=1,
        seed=1,
    )

    assert report.accepted_count == 1
    assert report.rejected_count == 5
    assert report.step_rows[0].accepted is True
    assert report.step_rows[0].proposal_valid is True
    assert report.step_rows[0].log_acceptance_ratio == 5.0
    assert all(step_row.proposal_valid is True for step_row in report.step_rows)
    assert any(not step_row.accepted for step_row in report.step_rows)
    assert (
        rooted_topology_fingerprint(report.final_state.tree.to_tree())
        == _preferred_nni_topology_fingerprint()
    )


def test_nni_topology_proposal_rejects_non_bifurcating_tree() -> None:
    current_state = score_bayesian_phylogenetic_state(
        tree=PhyloTree(
            TreeNode(
                children=[
                    TreeNode(name="A", branch_length=0.4),
                    TreeNode(name="B", branch_length=0.5),
                    TreeNode(name="C", branch_length=0.6),
                ]
            ),
            rooted=True,
        ),
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    proposal = propose_nni_topology_move(current_state, Random(7))

    assert proposal.is_valid is False
    assert proposal.invalid_reason == ("rooted NNI enumeration requires a binary root")
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def _build_scored_nni_state(
    *,
    update_prior_components,
) -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
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
        update_prior_components=update_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )


def _preferred_nni_topology_fingerprint() -> str:
    preferred_tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(name="A", branch_length=0.4),
                TreeNode(
                    children=[
                        TreeNode(name="B", branch_length=0.5),
                        TreeNode(name="C", branch_length=1.0),
                    ],
                    branch_length=0.6,
                ),
            ]
        ),
        rooted=True,
    )
    return rooted_topology_fingerprint(preferred_tree)


def _topology_preference_prior_components(
    state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    topology_fingerprint = rooted_topology_fingerprint(state.tree.to_tree())
    return [
        BayesianPriorComponentState(
            component_name="topology-preference",
            family="discrete",
            log_prior=0.0
            if topology_fingerprint == _preferred_nni_topology_fingerprint()
            else -5.0,
        )
    ]


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
