from __future__ import annotations

import math
from pathlib import Path
from random import Random

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_reversible_jump_model_switch_move,
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
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
)
from bijux_phylogenetics.phylo.topology.clades import rooted_topology_fingerprint
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_reversible_jump_model_switch_records_jc69_to_k80_mapping_and_acceptance_ratio() -> (
    None
):
    current_state = _build_scored_jc69_state()

    proposal = propose_reversible_jump_model_switch_move(
        current_state,
        Random(2),
        log_kappa_standard_deviation=0.6,
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.changed_fields == (
        "categorical_parameters.substitution-model",
        "scalar_parameters.kappa",
    )
    assert proposal.proposed_tree is not None
    assert proposal.proposed_model_parameters is not None
    assert rooted_topology_fingerprint(
        proposal.proposed_tree
    ) == rooted_topology_fingerprint(current_state.tree.to_tree())
    assert (
        current_state.model_parameters.categorical_parameters["substitution-model"]
        == "JC69"
    )
    assert (
        proposal.proposed_model_parameters.categorical_parameters["substitution-model"]
        == "K80"
    )
    proposed_kappa = proposal.proposed_model_parameters.scalar_parameters["kappa"]
    assert proposed_kappa > 0.0
    assert "kappa" not in current_state.model_parameters.scalar_parameters
    assert math.isclose(
        proposal.log_forward_density,
        _lognormal_positive_draw_density(
            proposed_value=proposed_kappa,
            log_standard_deviation=0.6,
        ),
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert proposal.log_reverse_density == 0.0
    assert math.isclose(
        proposal.log_hastings_ratio,
        -proposal.log_forward_density,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_reversible_jump_model_switch_records_k80_to_jc69_mapping_and_acceptance_ratio() -> (
    None
):
    current_state = _build_scored_k80_state(kappa=2.5)

    proposal = propose_reversible_jump_model_switch_move(
        current_state,
        Random(5),
        log_kappa_standard_deviation=0.6,
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.changed_fields == (
        "categorical_parameters.substitution-model",
        "scalar_parameters.kappa",
    )
    assert proposal.proposed_model_parameters is not None
    assert (
        current_state.model_parameters.categorical_parameters["substitution-model"]
        == "K80"
    )
    assert (
        proposal.proposed_model_parameters.categorical_parameters["substitution-model"]
        == "JC69"
    )
    assert "kappa" not in proposal.proposed_model_parameters.scalar_parameters
    assert proposal.log_forward_density == 0.0
    assert math.isclose(
        proposal.log_reverse_density,
        _lognormal_positive_draw_density(
            proposed_value=2.5,
            log_standard_deviation=0.6,
        ),
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        proposal.log_hastings_ratio,
        proposal.log_reverse_density,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_reversible_jump_model_switch_changes_real_sequence_likelihood_surface() -> (
    None
):
    current_state = _build_scored_jc69_state()

    proposal = propose_reversible_jump_model_switch_move(
        current_state,
        Random(3),
        log_kappa_standard_deviation=0.45,
    )

    assert proposal.is_valid is True
    assert proposal.proposed_tree is not None
    assert proposal.proposed_model_parameters is not None

    proposed_state = score_bayesian_phylogenetic_state(
        tree=proposal.proposed_tree,
        model_parameters=proposal.proposed_model_parameters,
        update_prior_components=_flat_prior_components,
        update_log_likelihood=_sequence_model_log_likelihood,
    )

    assert not math.isclose(
        proposed_state.log_likelihood,
        current_state.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_sampler_uses_reversible_jump_model_switch_on_real_sequence_surface() -> None:
    initial_state = _build_scored_jc69_state()

    report = run_metropolis_hastings_sampler(
        initial_state=initial_state,
        propose_state=lambda current_state, rng: (
            propose_reversible_jump_model_switch_move(
                current_state,
                rng,
                log_kappa_standard_deviation=0.45,
            )
        ),
        update_prior_components=_flat_prior_components,
        update_log_likelihood=_sequence_model_log_likelihood,
        iteration_count=8,
        sample_every=1,
        seed=11,
    )

    sampled_models = [
        state.model_parameters.categorical_parameters["substitution-model"]
        for state in report.sampled_states
    ]
    assert report.accepted_count >= 1
    assert all(step.proposal_valid is True for step in report.step_rows)
    assert all(
        step.proposal_changed_fields
        == (
            "categorical_parameters.substitution-model",
            "scalar_parameters.kappa",
        )
        for step in report.step_rows
    )
    assert "JC69" in sampled_models
    assert "K80" in sampled_models
    assert any(
        step.log_acceptance_ratio is not None
        and not math.isclose(
            step.log_acceptance_ratio,
            step.log_hastings_ratio,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for step in report.step_rows
    )


def test_reversible_jump_model_switch_requires_substitution_model_label() -> None:
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
        update_prior_components=_flat_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    proposal = propose_reversible_jump_model_switch_move(
        current_state,
        Random(0),
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "reversible-jump model-switch proposal requires one "
        "'substitution-model' categorical parameter"
    )
    assert proposal.proposed_tree is None
    assert proposal.proposed_model_parameters is None


def test_reversible_jump_model_switch_rejects_inconsistent_jc69_state_with_kappa() -> (
    None
):
    current_state = score_bayesian_phylogenetic_state(
        tree=load_tree(fixture("trees", "k80_likelihood_tree_2_taxa.nwk")),
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "JC69"},
            scalar_parameters={"kappa": 2.0},
        ),
        update_prior_components=_flat_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    proposal = propose_reversible_jump_model_switch_move(
        current_state,
        Random(0),
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "reversible-jump model-switch proposal requires JC69 states to omit "
        "the standalone 'kappa' scalar parameter"
    )


def test_reversible_jump_model_switch_rejects_inconsistent_k80_state_without_kappa() -> (
    None
):
    current_state = score_bayesian_phylogenetic_state(
        tree=load_tree(fixture("trees", "k80_likelihood_tree_2_taxa.nwk")),
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "K80"}
        ),
        update_prior_components=_flat_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    proposal = propose_reversible_jump_model_switch_move(
        current_state,
        Random(0),
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "reversible-jump model-switch proposal requires K80 states to include "
        "one positive 'kappa' scalar parameter"
    )


def test_reversible_jump_model_switch_rejects_unsupported_model_label() -> None:
    current_state = score_bayesian_phylogenetic_state(
        tree=load_tree(fixture("trees", "k80_likelihood_tree_2_taxa.nwk")),
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "HKY85"},
            scalar_parameters={"kappa": 2.0},
        ),
        update_prior_components=_flat_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    proposal = propose_reversible_jump_model_switch_move(
        current_state,
        Random(0),
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "reversible-jump model-switch proposal supports only JC69 and K80 "
        "within the nucleotide-substitution-model family"
    )


def _build_scored_jc69_state() -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=load_tree(fixture("trees", "k80_likelihood_tree_2_taxa.nwk")),
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "JC69"}
        ),
        update_prior_components=_flat_prior_components,
        update_log_likelihood=_sequence_model_log_likelihood,
    )


def _build_scored_k80_state(*, kappa: float) -> BayesianPhylogeneticState:
    return score_bayesian_phylogenetic_state(
        tree=load_tree(fixture("trees", "k80_likelihood_tree_2_taxa.nwk")),
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters={"substitution-model": "K80"},
            scalar_parameters={"kappa": kappa},
        ),
        update_prior_components=_flat_prior_components,
        update_log_likelihood=_sequence_model_log_likelihood,
    )


def _sequence_model_log_likelihood(state: BayesianPhylogeneticState) -> float:
    model_name = state.model_parameters.categorical_parameters["substitution-model"]
    tree = state.tree.to_tree()
    tree.rooted = state.tree.rooted
    records = load_fasta_alignment(
        fixture("alignments", "k80_likelihood_alignment_2_taxa.fasta")
    )
    if model_name == "JC69":
        return evaluate_jc69_tree_likelihood(tree, records).log_likelihood
    if model_name == "K80":
        return evaluate_k80_tree_likelihood(
            tree,
            records,
            kappa=state.model_parameters.scalar_parameters["kappa"],
        ).log_likelihood
    raise AssertionError(f"unexpected substitution model in test state: {model_name}")


def _flat_prior_components(
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


def _lognormal_positive_draw_density(
    *,
    proposed_value: float,
    log_standard_deviation: float,
) -> float:
    return (
        -math.log(proposed_value)
        - math.log(log_standard_deviation)
        - (0.5 * math.log(2.0 * math.pi))
        - (math.log(proposed_value) ** 2 / (2.0 * (log_standard_deviation**2)))
    )
