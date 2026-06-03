from __future__ import annotations

import math
from pathlib import Path
from random import Random

import numpy
import pytest

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_gamma_alpha_move,
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
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma,
)
from bijux_phylogenetics.phylo.likelihood.gamma import (
    build_discrete_gamma_rate_categories,
)
from bijux_phylogenetics.phylo.topology.clades import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_gamma_alpha_proposal_keeps_alpha_on_positive_support() -> None:
    current_state = _build_scored_gamma_alpha_state()
    rng = Random(17)

    for _ in range(100):
        proposal = propose_gamma_alpha_move(
            current_state,
            rng,
            log_scale_standard_deviation=0.6,
        )
        assert proposal.is_valid is True
        assert proposal.invalid_reason is None
        assert proposal.changed_fields == ("scalar_parameters.gamma-alpha",)
        assert proposal.proposed_model_parameters is not None
        assert proposal.proposed_model_parameters.scalar_parameters["gamma-alpha"] > 0.0


def test_gamma_alpha_proposal_changes_category_rates_and_likelihood() -> None:
    current_state = _build_scored_gamma_alpha_state()

    proposal = propose_gamma_alpha_move(
        current_state,
        Random(2),
        log_scale_standard_deviation=0.6,
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.changed_fields == ("scalar_parameters.gamma-alpha",)
    proposed_tree = proposal.proposed_tree
    proposed_model_parameters = proposal.proposed_model_parameters
    assert proposed_tree is not None
    assert proposed_model_parameters is not None
    assert rooted_topology_fingerprint(proposed_tree) == rooted_topology_fingerprint(
        current_state.tree.to_tree()
    )
    current_alpha = current_state.model_parameters.scalar_parameters["gamma-alpha"]
    proposed_alpha = proposed_model_parameters.scalar_parameters["gamma-alpha"]
    assert proposed_alpha > 0.0
    assert math.isclose(
        proposal.log_hastings_ratio,
        math.log(proposed_alpha / current_alpha),
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert [
        round(category.rate, 12)
        for category in build_discrete_gamma_rate_categories(
            alpha=current_alpha,
            category_count=4,
        )
    ] != [
        round(category.rate, 12)
        for category in build_discrete_gamma_rate_categories(
            alpha=proposed_alpha,
            category_count=4,
        )
    ]

    proposed_state = score_bayesian_phylogenetic_state(
        tree=proposed_tree,
        model_parameters=proposed_model_parameters,
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_gamma_log_likelihood,
    )

    assert not math.isclose(
        proposed_state.log_likelihood,
        current_state.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


@pytest.mark.slow
def test_sampler_uses_gamma_alpha_proposal_on_real_heterogeneity_surface() -> None:
    initial_state = _build_scored_gamma_alpha_state()

    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: propose_gamma_alpha_move(
            current_state,
            rng,
            log_scale_standard_deviation=0.6,
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_gamma_log_likelihood,
        iteration_count=6,
        sample_every=1,
        seed=0,
    )

    assert report.accepted_count == 5
    assert report.rejected_count == 1
    assert all(step_row.proposal_valid is True for step_row in report.step_rows)
    assert all(
        step_row.proposal_changed_fields == ("scalar_parameters.gamma-alpha",)
        for step_row in report.step_rows
    )
    assert any(
        step_row.log_acceptance_ratio is not None
        and not math.isclose(
            step_row.log_acceptance_ratio,
            step_row.log_hastings_ratio,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for step_row in report.step_rows
    )
    assert report.final_state.model_parameters.scalar_parameters["gamma-alpha"] > 0.0


def test_gamma_alpha_proposal_requires_gamma_alpha_scalar() -> None:
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

    proposal = propose_gamma_alpha_move(
        current_state,
        Random(0),
        log_scale_standard_deviation=0.5,
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "gamma-alpha proposal requires one 'gamma-alpha' scalar parameter"
    )
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def _build_scored_gamma_alpha_state() -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=load_tree(
            fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk")
        ),
        model_parameters=build_bayesian_model_parameter_state(
            scalar_parameters={"gamma-alpha": 0.8}
        ),
        update_prior_components=_zero_prior_components,
        update_log_likelihood=_gamma_log_likelihood,
    )


def _gamma_log_likelihood(state: BayesianPhylogeneticState) -> float:
    return evaluate_empirical_protein_tree_likelihood_with_discrete_gamma(
        load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk")),
        load_fasta_alignment(
            fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta")
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        alpha=state.model_parameters.scalar_parameters["gamma-alpha"],
        category_count=4,
        matrix_label="compact-polar",
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


def _compact_polar_rate_matrix() -> numpy.ndarray:
    state_order = _protein_state_order()
    state_index = {state: index for index, state in enumerate(state_order)}
    rate_matrix = numpy.full((len(state_order), len(state_order)), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for (left_state, right_state), rate in {
        ("A", "C"): 0.45,
        ("C", "D"): 0.35,
        ("D", "E"): 0.55,
        ("A", "E"): 0.20,
    }.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(rate_matrix[row_index, :].sum())
    return rate_matrix


def _biased_root_prior() -> numpy.ndarray:
    prior = numpy.full(20, 0.02, dtype=float)
    state_index = {state: index for index, state in enumerate(_protein_state_order())}
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())


def _protein_state_order() -> tuple[str, ...]:
    return (
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    )
