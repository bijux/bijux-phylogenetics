from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.prior_sampling import (
    sample_prior_only_phylogenetic_state,
)
from bijux_phylogenetics.bayesian.state import (
    build_bayesian_model_parameter_state,
    build_bayesian_phylogenetic_state,
    build_bayesian_phylogenetic_state_from_prior_only_sample,
    build_bayesian_prior_component_state,
    deserialize_bayesian_phylogenetic_state,
    deserialize_bayesian_phylogenetic_state_json,
    serialize_bayesian_phylogenetic_state,
    serialize_bayesian_phylogenetic_state_json,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_beta_probability_substitution_parameter_prior,
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    build_uniform_rooted_tree_topology_prior,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    DNA_EXCHANGEABILITY_LABELS,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_bayesian_state_round_trip_preserves_posterior_score_exactly() -> None:
    sample = sample_prior_only_phylogenetic_state(
        tree_topology_prior=build_uniform_rooted_tree_topology_prior(
            ["A", "B", "C", "D"]
        ),
        branch_length_prior=build_exponential_branch_length_prior(rate=1.25),
        substitution_parameter_prior=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=0.5
            ),
            exchangeability_prior=build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=DNA_EXCHANGEABILITY_LABELS,
                concentration_parameters={
                    "AC": 1.5,
                    "AG": 3.0,
                    "AT": 2.0,
                    "CG": 1.0,
                    "CT": 2.5,
                    "GT": 1.25,
                },
            ),
            base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=("A", "C", "G", "T"),
                concentration_parameters={
                    "A": 3.0,
                    "C": 2.0,
                    "G": 1.5,
                    "T": 1.0,
                },
            ),
            invariant_proportion_prior=build_beta_probability_substitution_parameter_prior(
                alpha=2.0,
                beta=5.0,
            ),
        ),
        seed=19,
    )
    state = build_bayesian_phylogenetic_state_from_prior_only_sample(
        sample,
        log_likelihood=-12.375,
    )

    serialized_json = serialize_bayesian_phylogenetic_state_json(state)
    restored_state = deserialize_bayesian_phylogenetic_state_json(serialized_json)

    assert restored_state == state
    assert restored_state.posterior_log_score == state.posterior_log_score
    assert (
        restored_state.posterior_log_score
        == restored_state.total_log_prior + restored_state.log_likelihood
    )


def test_bayesian_state_deserialization_rejects_tampered_posterior_score() -> None:
    state = build_bayesian_phylogenetic_state(
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
            categorical_parameters={"substitution-model": "K80"},
            scalar_parameters={"kappa": 2.0},
        ),
        prior_components=[
            build_bayesian_prior_component_state(
                component_name="tree-topology",
                family="uniform-rooted-labeled-bifurcating",
                log_prior=-1.0,
            ),
            build_bayesian_prior_component_state(
                component_name="branch-lengths",
                family="exponential",
                log_prior=-0.7,
            ),
        ],
        log_likelihood=-4.25,
    )
    tampered_payload = serialize_bayesian_phylogenetic_state(state)
    tampered_payload["posterior_log_score"] = -999.0

    with pytest.raises(
        PhylogeneticsError,
        match="posterior_log_score does not match serialized payload",
    ):
        deserialize_bayesian_phylogenetic_state(tampered_payload)


def test_bayesian_state_requires_nonempty_prior_components() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="requires at least one prior component",
    ):
        build_bayesian_phylogenetic_state(
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
            prior_components=[],
            log_likelihood=0.0,
        )
