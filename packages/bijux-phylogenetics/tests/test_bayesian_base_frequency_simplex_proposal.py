from __future__ import annotations

import math
from pathlib import Path
from random import Random

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_base_frequency_simplex_move,
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
from bijux_phylogenetics.phylo.likelihood import evaluate_f81_tree_likelihood
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    parameterize_dna_base_frequency_simplex,
)
from bijux_phylogenetics.phylo.topology.clades import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_base_frequency_simplex_proposal_changes_likelihood_and_keeps_simplex_valid() -> (
    None
):
    current_state = _build_scored_base_frequency_state()

    proposal = propose_base_frequency_simplex_move(
        current_state,
        Random(2),
        unconstrained_coordinate_standard_deviation=0.6,
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.changed_fields == ("vector_parameters.base-frequencies.A",)
    assert math.isclose(
        proposal.log_hastings_ratio,
        0.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    proposed_tree = proposal.proposed_tree
    proposed_model_parameters = proposal.proposed_model_parameters
    assert proposed_tree is not None
    assert proposed_model_parameters is not None
    assert rooted_topology_fingerprint(proposed_tree) == rooted_topology_fingerprint(
        current_state.tree.to_tree()
    )
    assert (
        proposed_model_parameters.categorical_parameters["substitution-model"] == "F81"
    )
    proposed_base_frequencies = proposed_model_parameters.vector_parameters[
        "base-frequencies"
    ]
    assert tuple(sorted(proposed_base_frequencies)) == ("A", "C", "G", "T")
    assert all(value > 0.0 for value in proposed_base_frequencies.values())
    assert math.isclose(
        sum(proposed_base_frequencies.values()),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    proposed_state = score_bayesian_phylogenetic_state(
        tree=proposed_tree,
        model_parameters=proposed_model_parameters,
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_f81_log_likelihood,
    )

    assert not math.isclose(
        proposed_state.log_likelihood,
        current_state.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_sampler_uses_base_frequency_simplex_proposal_on_real_likelihood_surface() -> (
    None
):
    initial_state = _build_scored_base_frequency_state()

    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: propose_base_frequency_simplex_move(
            current_state,
            rng,
            unconstrained_coordinate_standard_deviation=0.6,
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_f81_log_likelihood,
        iteration_count=6,
        sample_every=1,
        seed=0,
    )

    assert report.accepted_count == 4
    assert report.rejected_count == 2
    assert all(step_row.proposal_valid is True for step_row in report.step_rows)
    assert all(
        changed_field.startswith("vector_parameters.base-frequencies.")
        for step_row in report.step_rows
        for changed_field in step_row.proposal_changed_fields
    )
    assert any(
        step_row.log_acceptance_ratio is not None
        and not math.isclose(
            step_row.log_acceptance_ratio,
            0.0,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for step_row in report.step_rows
    )
    final_base_frequencies = report.final_state.model_parameters.vector_parameters[
        "base-frequencies"
    ]
    assert all(value > 0.0 for value in final_base_frequencies.values())
    assert math.isclose(
        sum(final_base_frequencies.values()),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_base_frequency_simplex_proposal_requires_base_frequency_vector() -> None:
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
        model_parameters=build_bayesian_model_parameter_state(),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    proposal = propose_base_frequency_simplex_move(
        current_state,
        Random(0),
        unconstrained_coordinate_standard_deviation=0.5,
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "base-frequency simplex proposal requires one 'base-frequencies' vector parameter"
    )
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def _build_scored_base_frequency_state() -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=load_tree(fixture("trees", "f81_likelihood_tree_2_taxa.nwk")),
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "F81"},
            vector_parameters={
                "base-frequencies": parameterize_dna_base_frequency_simplex(
                    {
                        "A": 0.4,
                        "C": 0.1,
                        "G": 0.2,
                        "T": 0.3,
                    }
                ).constrained_mapping()
            },
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_f81_log_likelihood,
    )


def _f81_log_likelihood(state: BayesianPhylogeneticState) -> float:
    return evaluate_f81_tree_likelihood(
        state.tree.to_tree(),
        load_fasta_alignment(
            fixture("alignments", "f81_likelihood_alignment_2_taxa.fasta")
        ),
        base_frequencies=state.model_parameters.vector_parameters["base-frequencies"],
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
