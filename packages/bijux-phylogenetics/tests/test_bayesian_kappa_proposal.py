from __future__ import annotations

import math
from pathlib import Path
from random import Random

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_kappa_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import evaluate_k80_tree_likelihood
from bijux_phylogenetics.phylo.topology.clades import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_kappa_proposal_changes_real_likelihood_and_preserves_model_label() -> None:
    current_state = _build_scored_k80_state()

    proposal = propose_kappa_move(
        current_state,
        Random(2),
        log_scale_standard_deviation=0.6,
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.changed_fields == ("scalar_parameters.kappa",)
    proposed_tree = proposal.proposed_tree
    proposed_model_parameters = proposal.proposed_model_parameters
    assert proposed_tree is not None
    assert proposed_model_parameters is not None
    assert rooted_topology_fingerprint(proposed_tree) == rooted_topology_fingerprint(
        current_state.tree.to_tree()
    )
    assert (
        proposed_model_parameters.categorical_parameters["substitution-model"] == "K80"
    )
    current_kappa = current_state.model_parameters.scalar_parameters["kappa"]
    proposed_kappa = proposed_model_parameters.scalar_parameters["kappa"]
    assert proposed_kappa > 0.0
    assert math.isclose(
        proposal.log_hastings_ratio,
        math.log(proposed_kappa / current_kappa),
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    proposed_state = score_bayesian_phylogenetic_state(
        tree=proposed_tree,
        model_parameters=proposed_model_parameters,
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_k80_log_likelihood,
    )

    assert not math.isclose(
        proposed_state.log_likelihood,
        current_state.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_sampler_uses_kappa_proposal_on_real_k80_surface() -> None:
    initial_state = _build_scored_k80_state()

    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: propose_kappa_move(
            current_state,
            rng,
            log_scale_standard_deviation=0.6,
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_k80_log_likelihood,
        iteration_count=6,
        sample_every=1,
        seed=0,
    )

    assert report.accepted_count >= 1
    assert all(step_row.proposal_valid is True for step_row in report.step_rows)
    assert all(
        step_row.proposal_changed_fields == ("scalar_parameters.kappa",)
        for step_row in report.step_rows
    )
    assert all(
        state.model_parameters.categorical_parameters["substitution-model"] == "K80"
        for state in report.sampled_states
    )


def test_kappa_proposal_requires_scalar_parameter() -> None:
    current_state = score_bayesian_phylogenetic_state(
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
            categorical_parameters={"substitution-model": "K80"}
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    proposal = propose_kappa_move(
        current_state,
        Random(0),
        log_scale_standard_deviation=0.5,
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "kappa proposal requires one 'kappa' scalar parameter"
    )
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def _build_scored_k80_state() -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=load_tree(fixture("trees", "k80_likelihood_tree_2_taxa.nwk")),
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "K80"},
            scalar_parameters={"kappa": 2.0},
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_k80_log_likelihood,
    )


def _k80_log_likelihood(state: BayesianPhylogeneticState) -> float:
    return evaluate_k80_tree_likelihood(
        state.tree.to_tree(),
        load_fasta_alignment(
            fixture("alignments", "k80_likelihood_alignment_2_taxa.fasta")
        ),
        kappa=state.model_parameters.scalar_parameters["kappa"],
    ).log_likelihood


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
