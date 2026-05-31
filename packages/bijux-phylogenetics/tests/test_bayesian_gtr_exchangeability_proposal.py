from __future__ import annotations

import math
from pathlib import Path
from random import Random

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_gtr_exchangeability_move,
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
from bijux_phylogenetics.phylo.likelihood import evaluate_gtr_tree_likelihood
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    DNA_EXCHANGEABILITY_LABELS,
    parameterize_dna_exchangeability_simplex,
)
from bijux_phylogenetics.phylo.topology.clades import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

FIXTURES = Path(__file__).parent / "fixtures"
_GTR_BASE_FREQUENCIES = [0.4, 0.1, 0.2, 0.3]


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_gtr_exchangeability_proposal_changes_likelihood_and_keeps_simplex_valid() -> (
    None
):
    current_state = _build_scored_gtr_exchangeability_state()

    proposal = propose_gtr_exchangeability_move(
        current_state,
        Random(2),
        unconstrained_coordinate_standard_deviation=0.6,
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.changed_fields == ("vector_parameters.exchangeabilities.AC",)
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
        proposed_model_parameters.categorical_parameters["substitution-model"] == "GTR"
    )
    proposed_exchangeabilities = proposed_model_parameters.vector_parameters[
        "exchangeabilities"
    ]
    assert tuple(sorted(proposed_exchangeabilities)) == tuple(
        sorted(DNA_EXCHANGEABILITY_LABELS)
    )
    assert all(value > 0.0 for value in proposed_exchangeabilities.values())
    assert math.isclose(
        sum(proposed_exchangeabilities.values()),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    proposed_state = score_bayesian_phylogenetic_state(
        tree=proposed_tree,
        model_parameters=proposed_model_parameters,
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_gtr_log_likelihood,
    )

    assert not math.isclose(
        proposed_state.log_likelihood,
        current_state.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_sampler_uses_gtr_exchangeability_proposal_on_real_likelihood_surface() -> None:
    initial_state = _build_scored_gtr_exchangeability_state()

    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: propose_gtr_exchangeability_move(
            current_state,
            rng,
            unconstrained_coordinate_standard_deviation=0.6,
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_gtr_log_likelihood,
        iteration_count=6,
        sample_every=1,
        seed=0,
    )

    assert report.accepted_count == 5
    assert report.rejected_count == 1
    assert all(step_row.proposal_valid is True for step_row in report.step_rows)
    assert all(
        changed_field.startswith("vector_parameters.exchangeabilities.")
        for step_row in report.step_rows
        for changed_field in step_row.proposal_changed_fields
    )
    assert any(
        step_row.log_acceptance_ratio is not None
        and not math.isclose(
            step_row.log_acceptance_ratio, 0.0, rel_tol=0.0, abs_tol=1e-12
        )
        for step_row in report.step_rows
    )
    final_exchangeabilities = report.final_state.model_parameters.vector_parameters[
        "exchangeabilities"
    ]
    assert all(value > 0.0 for value in final_exchangeabilities.values())
    assert math.isclose(
        sum(final_exchangeabilities.values()),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_gtr_exchangeability_proposal_requires_exchangeability_vector() -> None:
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

    proposal = propose_gtr_exchangeability_move(
        current_state,
        Random(0),
        unconstrained_coordinate_standard_deviation=0.5,
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "GTR exchangeability proposal requires one 'exchangeabilities' vector parameter"
    )
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def _build_scored_gtr_exchangeability_state() -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=load_tree(fixture("trees", "gtr_likelihood_tree_2_taxa.nwk")),
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "GTR"},
            vector_parameters={
                "exchangeabilities": parameterize_dna_exchangeability_simplex(
                    {
                        "AC": 1.0,
                        "AG": 2.0,
                        "AT": 0.5,
                        "CG": 1.5,
                        "CT": 1.75,
                        "GT": 1.25,
                    }
                ).constrained_mapping()
            },
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_gtr_log_likelihood,
    )


def _gtr_log_likelihood(state: BayesianPhylogeneticState) -> float:
    return evaluate_gtr_tree_likelihood(
        state.tree.to_tree(),
        load_fasta_alignment(
            fixture("alignments", "gtr_likelihood_alignment_2_taxa.fasta")
        ),
        exchangeabilities=state.model_parameters.vector_parameters["exchangeabilities"],
        base_frequencies=_GTR_BASE_FREQUENCIES,
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
