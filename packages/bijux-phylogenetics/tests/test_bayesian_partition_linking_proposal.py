from __future__ import annotations

import math
from random import Random

from bijux_phylogenetics.bayesian.metropolis_hastings import (
    propose_partition_linking_move,
    run_metropolis_hastings_sampler,
    score_bayesian_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.partition_model_priors import (
    PartitionModelPriorEvaluationReport,
    PartitionSubstitutionParameterState,
    build_partition_model_prior_bundle,
    build_partition_parameter_linkage_plan,
    build_partition_substitution_model_definition,
    evaluate_partition_model_log_prior,
)
from bijux_phylogenetics.bayesian.partition_model_state import (
    build_partition_model_parameter_state,
    resolve_partition_parameter_linkage_plan_from_model_parameters,
    resolve_partition_parameter_states_from_model_parameters,
)
from bijux_phylogenetics.bayesian.state import (
    BayesianModelParameterState,
    BayesianPhylogeneticState,
    BayesianPriorComponentState,
    build_bayesian_model_parameter_state,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def test_partition_linking_move_unlinks_one_shared_kappa_parameter() -> None:
    current_state = _build_scored_partition_state(
        kappa_linkage_policy="linked",
        kappa_values=(2.0, 2.0),
    )
    current_report = _partition_prior_report(current_state.model_parameters)

    proposal = propose_partition_linking_move(
        current_state,
        Random(0),
        partition_models=_partition_models(),
        target_name="kappa",
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.proposed_model_parameters is not None

    proposed_report = _partition_prior_report(proposal.proposed_model_parameters)
    sampler_report = run_metropolis_hastings_sampler(
        initial_state=current_state,
        propose_state=lambda state, rng: propose_partition_linking_move(
            state,
            rng,
            partition_models=_partition_models(),
            target_name="kappa",
        ),
        update_prior_components=_partition_prior_components,
        update_log_likelihood=_zero_log_likelihood,
        iteration_count=1,
        sample_every=1,
        seed=13,
    )

    assert current_report.parameter_count == 3
    assert proposed_report.parameter_count == 4
    assert (
        proposal.proposed_model_parameters.categorical_parameters[
            "partition-linkage:kappa:gene_alpha"
        ]
        == "gene_alpha"
    )
    assert (
        proposal.proposed_model_parameters.categorical_parameters[
            "partition-linkage:kappa:gene_beta"
        ]
        == "gene_beta"
    )
    assert "partition-parameter:kappa:gene_alpha" in (
        proposal.proposed_model_parameters.scalar_parameters
    )
    assert "partition-parameter:kappa:gene_beta" in (
        proposal.proposed_model_parameters.scalar_parameters
    )
    assert "partition-parameter:kappa:kappa-shared" not in (
        proposal.proposed_model_parameters.scalar_parameters
    )
    assert not math.isclose(
        sampler_report.step_rows[0].proposed_posterior_log_score or 0.0,
        sampler_report.step_rows[0].current_posterior_log_score,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_partition_linking_move_links_matching_unlinked_kappa_parameters() -> None:
    current_state = _build_scored_partition_state(
        kappa_linkage_policy="unlinked",
        kappa_values=(2.0, 2.0),
    )
    current_report = _partition_prior_report(current_state.model_parameters)

    proposal = propose_partition_linking_move(
        current_state,
        Random(0),
        partition_models=_partition_models(),
        target_name="kappa",
    )

    assert proposal.is_valid is True
    assert proposal.invalid_reason is None
    assert proposal.proposed_model_parameters is not None

    proposed_report = _partition_prior_report(proposal.proposed_model_parameters)
    proposed_state = score_bayesian_phylogenetic_state(
        tree=proposal.proposed_tree,
        model_parameters=proposal.proposed_model_parameters,
        update_prior_components=_partition_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )

    assert current_report.parameter_count == 4
    assert proposed_report.parameter_count == 3
    assert (
        proposal.proposed_model_parameters.categorical_parameters[
            "partition-linkage:kappa:gene_alpha"
        ]
        == "kappa-shared"
    )
    assert (
        proposal.proposed_model_parameters.categorical_parameters[
            "partition-linkage:kappa:gene_beta"
        ]
        == "kappa-shared"
    )
    assert (
        proposal.proposed_model_parameters.scalar_parameters[
            "partition-parameter:kappa:kappa-shared"
        ]
        == 2.0
    )
    assert not math.isclose(
        proposed_state.posterior_log_score,
        current_state.posterior_log_score,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_partition_linking_move_rejects_linking_mismatched_partition_values() -> None:
    current_state = _build_scored_partition_state(
        kappa_linkage_policy="unlinked",
        kappa_values=(2.0, 3.0),
    )

    proposal = propose_partition_linking_move(
        current_state,
        Random(0),
        partition_models=_partition_models(),
        target_name="kappa",
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "partition model state encoding requires linked groups to share one realized value"
    )


def test_partition_linking_move_rejects_mixed_target_linkage_states() -> None:
    partition_models = (
        build_partition_substitution_model_definition(
            partition_name="gene_alpha",
            model_name="HKY85+G",
        ),
        build_partition_substitution_model_definition(
            partition_name="gene_beta",
            model_name="HKY85+G",
        ),
        build_partition_substitution_model_definition(
            partition_name="gene_gamma",
            model_name="HKY85+G",
        ),
    )
    linked_model_parameters = build_partition_model_parameter_state(
        partition_models=partition_models,
        linkage_plan=build_partition_parameter_linkage_plan(
            partition_names=("gene_alpha", "gene_beta", "gene_gamma"),
            linkage_policies={
                "kappa": "linked",
                "base-frequencies": "linked",
                "gamma-alpha": "linked",
            },
        ),
        partition_parameter_states=(
            PartitionSubstitutionParameterState(
                partition_name="gene_alpha",
                kappa=2.0,
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
                gamma_alpha=0.7,
            ),
            PartitionSubstitutionParameterState(
                partition_name="gene_beta",
                kappa=2.0,
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
                gamma_alpha=0.7,
            ),
            PartitionSubstitutionParameterState(
                partition_name="gene_gamma",
                kappa=2.0,
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
                gamma_alpha=0.7,
            ),
        ),
        preserved_categorical_parameters={"substitution-model": "partitioned-dna"},
    )
    categorical_parameters = dict(linked_model_parameters.categorical_parameters)
    scalar_parameters = dict(linked_model_parameters.scalar_parameters)
    categorical_parameters["partition-linkage:kappa:gene_gamma"] = "gene_gamma"
    scalar_parameters["partition-parameter:kappa:gene_gamma"] = 2.0
    current_state = score_bayesian_phylogenetic_state(
        tree=_test_tree(),
        model_parameters=build_bayesian_model_parameter_state(
            categorical_parameters=categorical_parameters,
            scalar_parameters=scalar_parameters,
            vector_parameters=linked_model_parameters.vector_parameters,
        ),
        update_prior_components=lambda state: _partition_prior_components_for_models(
            state,
            partition_models,
        ),
        update_log_likelihood=_zero_log_likelihood,
    )

    proposal = propose_partition_linking_move(
        current_state,
        Random(0),
        partition_models=partition_models,
        target_name="kappa",
    )

    assert proposal.is_valid is False
    assert proposal.invalid_reason == (
        "partition-linking proposal currently supports only fully linked or fully "
        "unlinked states for one target"
    )


def _build_scored_partition_state(
    *,
    kappa_linkage_policy: str,
    kappa_values: tuple[float, float],
) -> BayesianPhylogeneticState:
    partition_models = _partition_models()
    linkage_plan = build_partition_parameter_linkage_plan(
        partition_names=("gene_alpha", "gene_beta"),
        linkage_policies={
            "kappa": kappa_linkage_policy,
            "base-frequencies": "linked",
            "gamma-alpha": "linked",
        },
    )
    model_parameters = build_partition_model_parameter_state(
        partition_models=partition_models,
        linkage_plan=linkage_plan,
        partition_parameter_states=(
            PartitionSubstitutionParameterState(
                partition_name="gene_alpha",
                kappa=kappa_values[0],
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
                gamma_alpha=0.7,
            ),
            PartitionSubstitutionParameterState(
                partition_name="gene_beta",
                kappa=kappa_values[1],
                base_frequencies=(0.3, 0.2, 0.1, 0.4),
                gamma_alpha=0.7,
            ),
        ),
        preserved_categorical_parameters={"substitution-model": "partitioned-dna"},
    )
    return score_bayesian_phylogenetic_state(
        tree=_test_tree(),
        model_parameters=model_parameters,
        update_prior_components=_partition_prior_components,
        update_log_likelihood=_zero_log_likelihood,
    )


def _partition_models():
    return (
        build_partition_substitution_model_definition(
            partition_name="gene_alpha",
            model_name="HKY85+G",
        ),
        build_partition_substitution_model_definition(
            partition_name="gene_beta",
            model_name="HKY85+G",
        ),
    )


def _partition_prior_report(
    model_parameters: BayesianModelParameterState,
) -> PartitionModelPriorEvaluationReport:
    return _partition_prior_report_for_models(model_parameters, _partition_models())


def _partition_prior_report_for_models(
    model_parameters: BayesianModelParameterState,
    partition_models,
) -> PartitionModelPriorEvaluationReport:
    linkage_plan = resolve_partition_parameter_linkage_plan_from_model_parameters(
        model_parameters=model_parameters,
        partition_names=tuple(
            partition_model.partition_name for partition_model in partition_models
        ),
    )
    partition_parameter_states = (
        resolve_partition_parameter_states_from_model_parameters(
            model_parameters=model_parameters,
            partition_models=partition_models,
            linkage_plan=linkage_plan,
        )
    )
    return evaluate_partition_model_log_prior(
        prior_bundle=build_partition_model_prior_bundle(
            partition_models=partition_models,
            linkage_plan=linkage_plan,
            substitution_prior_bundle=build_substitution_parameter_prior_bundle(
                kappa_prior=build_exponential_positive_substitution_parameter_prior(
                    rate=0.4
                ),
                base_frequency_prior=(
                    build_dirichlet_simplex_substitution_parameter_prior(
                        expected_component_names=("A", "C", "G", "T"),
                        concentration_parameters=(3.0, 2.0, 2.0, 3.0),
                    )
                ),
                gamma_alpha_prior=(
                    build_exponential_positive_substitution_parameter_prior(rate=0.9)
                ),
            ),
        ),
        partition_parameter_states=partition_parameter_states,
    )


def _partition_prior_components(
    state: BayesianPhylogeneticState,
) -> list[BayesianPriorComponentState]:
    return _partition_prior_components_for_models(state, _partition_models())


def _partition_prior_components_for_models(
    state: BayesianPhylogeneticState,
    partition_models,
) -> list[BayesianPriorComponentState]:
    report = _partition_prior_report_for_models(
        state.model_parameters,
        partition_models,
    )
    return [
        BayesianPriorComponentState(
            component_name="partition-model-prior",
            family="partition-model",
            log_prior=report.total_log_prior,
            parameter_values={"parameter-count": float(report.parameter_count)},
        )
    ]


def _test_tree() -> PhyloTree:
    return PhyloTree(
        TreeNode(
            children=[
                TreeNode(name="A", branch_length=0.1),
                TreeNode(name="B", branch_length=0.2),
            ]
        ),
        rooted=True,
    )


def _zero_log_likelihood(_state: BayesianPhylogeneticState) -> float:
    return 0.0
