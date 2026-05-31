from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_reversible_jump_model_switch_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.posterior_model_averaging import (
    summarize_metropolis_hastings_model_averaged_estimates,
    summarize_posterior_model_averaged_estimates,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
    build_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_jc69_tree_likelihood,
    evaluate_k80_tree_likelihood,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_posterior_model_averaging_emits_per_model_and_model_averaged_estimates() -> (
    None
):
    report = summarize_posterior_model_averaged_estimates(
        sampled_states=[
            _build_scored_jc69_state(),
            _build_scored_jc69_state(),
            _build_scored_k80_state(kappa=3.0),
        ]
    )

    assert report.sample_count == 3
    assert report.model_family == "nucleotide-substitution-model"
    assert report.sampled_models == ["JC69", "K80"]
    assert report.warnings == []

    support_by_model = {row.model_name: row for row in report.support_rows}
    assert math.isclose(
        support_by_model["JC69"].posterior_probability,
        2.0 / 3.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        support_by_model["K80"].posterior_probability,
        1.0 / 3.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    estimate_by_model = {row.model_name: row for row in report.per_model_estimate_rows}
    assert estimate_by_model["JC69"].estimate_name == "transition-transversion-ratio"
    assert estimate_by_model["JC69"].posterior_mean == 1.0
    assert estimate_by_model["JC69"].hpd_95_lower == 1.0
    assert estimate_by_model["JC69"].hpd_95_upper == 1.0
    assert estimate_by_model["K80"].posterior_mean == 3.0

    model_averaged_row = report.model_averaged_estimate_rows[0]
    assert model_averaged_row.estimate_name == "transition-transversion-ratio"
    assert math.isclose(
        model_averaged_row.posterior_mean,
        5.0 / 3.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert model_averaged_row.hpd_95_lower == 1.0
    assert model_averaged_row.hpd_95_upper == 3.0
    assert model_averaged_row.contributing_models == ["JC69", "K80"]


def test_metropolis_hastings_model_averaging_summarizes_real_reversible_jump_chain() -> (
    None
):
    initial_state = _build_scored_jc69_state()

    run_report = run_metropolis_hastings_sampler(
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

    report = summarize_metropolis_hastings_model_averaged_estimates(
        run_report=run_report
    )

    assert set(report.sampled_models) == {"JC69", "K80"}
    support_by_model = {row.model_name: row for row in report.support_rows}
    assert support_by_model["JC69"].supporting_sample_count >= 1
    assert support_by_model["K80"].supporting_sample_count >= 1
    estimate_by_model = {
        row.model_name: row.posterior_mean for row in report.per_model_estimate_rows
    }
    model_averaged_mean = report.model_averaged_estimate_rows[0].posterior_mean
    assert model_averaged_mean != estimate_by_model["JC69"]
    assert (
        min(estimate_by_model.values())
        <= model_averaged_mean
        <= max(estimate_by_model.values())
    )


def test_posterior_model_averaging_rejects_inconsistent_k80_state_without_kappa() -> (
    None
):
    with pytest.raises(PhylogeneticsError) as error_info:
        summarize_posterior_model_averaged_estimates(
            sampled_states=[
                build_bayesian_phylogenetic_state(
                    tree=load_tree(fixture("trees", "k80_likelihood_tree_2_taxa.nwk")),
                    model_parameters=build_bayesian_model_parameter_state(
                        categorical_parameters={"substitution-model": "K80"}
                    ),
                    prior_components=[
                        BayesianPriorComponentState(
                            component_name="flat-prior",
                            family="constant",
                            log_prior=0.0,
                        )
                    ],
                    log_likelihood=-1.0,
                )
            ]
        )

    assert error_info.value.code == "posterior_model_averaging_k80_kappa_missing"


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
