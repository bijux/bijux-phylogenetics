from __future__ import annotations

from random import Random

import pytest

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_tbr_topology_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.phylo.topology.clades import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.rooted_spr import enumerate_rooted_spr_neighbors
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def test_tbr_topology_proposal_reaches_tbr_only_neighbor_without_taxon_loss() -> None:
    current_state = _build_scored_tbr_state(
        update_prior_components=_zero_prior_components
    )
    current_tree = current_state.tree.to_tree()
    current_topology_fingerprint = rooted_topology_fingerprint(current_tree)
    spr_topology_fingerprints = {
        row.neighbor_topology_fingerprint
        for row in enumerate_rooted_spr_neighbors(current_tree).neighbor_rows
    }

    proposal = propose_tbr_topology_move(current_state, Random(2))

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.changed_fields == ("tree.topology:tbr:A|B|C|D:A:interface",)
    proposed_tree = proposal.proposed_tree
    assert proposed_tree is not None
    proposed_topology_fingerprint = rooted_topology_fingerprint(proposed_tree)
    assert set(proposed_tree.tip_names) == set(current_tree.tip_names)
    assert all(
        len(node.children) == 2
        for node in proposed_tree.iter_internal_nodes(order="preorder")
    )
    assert all(
        child.branch_length is not None for _parent, child in proposed_tree.iter_edges()
    )
    assert proposed_topology_fingerprint != current_topology_fingerprint
    assert proposed_topology_fingerprint not in spr_topology_fingerprints


@pytest.mark.slow
def test_tbr_topology_sampler_scores_and_accepts_rejects_moves() -> None:
    initial_state = _build_scored_tbr_state(
        update_prior_components=_topology_preference_prior_components
    )

    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=propose_tbr_topology_move,
        update_prior_components=_topology_preference_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=6,
        sample_every=1,
        seed=2,
    )

    assert report.accepted_count == 1
    assert report.rejected_count == 5
    assert report.step_rows[0].accepted is True
    assert report.step_rows[0].proposal_valid is True
    assert report.step_rows[0].log_acceptance_ratio == 5.0
    assert any(not step_row.accepted for step_row in report.step_rows)
    assert (
        rooted_topology_fingerprint(report.final_state.tree.to_tree())
        == _preferred_tbr_topology_fingerprint()
    )


def test_tbr_topology_proposal_rejects_non_bifurcating_tree() -> None:
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

    proposal = propose_tbr_topology_move(current_state, Random(7))

    assert proposal.is_valid is False
    assert proposal.invalid_reason == ("rooted TBR enumeration requires a binary root")
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def _build_scored_tbr_state(
    *,
    update_prior_components,
) -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=_build_tbr_start_tree(),
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=update_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )


def _build_tbr_start_tree() -> PhyloTree:
    tree = PhyloTree.from_newick("((((A:0.1,D:0.1):0.2,B:0.3):0.4,C:0.7):0.8,E:1.5);")
    tree.rooted = True
    return tree


def _preferred_tbr_topology_fingerprint() -> str:
    current_state = _build_scored_tbr_state(
        update_prior_components=_zero_prior_components
    )
    proposal = propose_tbr_topology_move(current_state, Random(2))
    proposed_tree = proposal.proposed_tree
    assert proposed_tree is not None
    return rooted_topology_fingerprint(proposed_tree)


def _topology_preference_prior_components(
    state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    topology_fingerprint = rooted_topology_fingerprint(state.tree.to_tree())
    return [
        BayesianPriorComponentState(
            component_name="topology-preference",
            family="discrete",
            log_prior=0.0
            if topology_fingerprint == _preferred_tbr_topology_fingerprint()
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
