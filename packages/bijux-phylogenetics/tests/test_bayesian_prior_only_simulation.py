from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
    build_fixed_branch_length_prior,
)
from bijux_phylogenetics.bayesian.prior_sampling import (
    PriorOnlySubstitutionParameterState,
    sample_prior_only_phylogenetic_state,
    simulate_prior_only_phylogenetic_states,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_beta_probability_substitution_parameter_prior,
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_fixed_positive_substitution_parameter_prior,
    build_fixed_simplex_substitution_parameter_prior,
    build_gamma_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    build_uniform_rooted_tree_topology_prior,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    DNA_EXCHANGEABILITY_LABELS,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_seeded_prior_only_simulation_is_reproducible_without_alignment() -> None:
    topology_prior = build_uniform_rooted_tree_topology_prior(["A", "B", "C", "D"])
    branch_length_prior = build_exponential_branch_length_prior(rate=1.5)
    substitution_prior = build_substitution_parameter_prior_bundle(
        kappa_prior=build_gamma_positive_substitution_parameter_prior(
            shape=2.5,
            scale=0.8,
        ),
        exchangeability_prior=build_dirichlet_simplex_substitution_parameter_prior(
            expected_component_names=DNA_EXCHANGEABILITY_LABELS,
            concentration_parameters={
                "AC": 2.0,
                "AG": 3.0,
                "AT": 2.0,
                "CG": 1.5,
                "CT": 2.5,
                "GT": 1.0,
            },
        ),
        base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
            expected_component_names=("A", "C", "G", "T"),
            concentration_parameters={
                "A": 4.0,
                "C": 2.0,
                "G": 3.0,
                "T": 1.5,
            },
        ),
        gamma_alpha_prior=build_exponential_positive_substitution_parameter_prior(
            rate=0.75
        ),
        invariant_proportion_prior=build_beta_probability_substitution_parameter_prior(
            alpha=2.0,
            beta=6.0,
        ),
    )

    left_report = simulate_prior_only_phylogenetic_states(
        tree_topology_prior=topology_prior,
        branch_length_prior=branch_length_prior,
        substitution_parameter_prior=substitution_prior,
        sample_count=2,
        seed=17,
    )
    right_report = simulate_prior_only_phylogenetic_states(
        tree_topology_prior=topology_prior,
        branch_length_prior=branch_length_prior,
        substitution_parameter_prior=substitution_prior,
        sample_count=2,
        seed=17,
    )
    changed_report = simulate_prior_only_phylogenetic_states(
        tree_topology_prior=topology_prior,
        branch_length_prior=branch_length_prior,
        substitution_parameter_prior=substitution_prior,
        sample_count=2,
        seed=23,
    )

    assert left_report == right_report
    assert changed_report != left_report
    assert left_report.sample_count == 2
    assert left_report.substitution_prior_count == 5
    assert len(left_report.samples) == 2
    assert left_report.samples[0].branch_rows
    assert left_report.samples[0].substitution_prior_rows
    assert left_report.samples[0].substitution_parameter_state.kappa is not None
    assert math.isclose(
        left_report.samples[0].total_log_prior,
        left_report.samples[0].topology_log_prior
        + left_report.samples[0].branch_length_log_prior
        + left_report.samples[0].substitution_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_prior_only_simulation_surfaces_fixed_parameter_samples() -> None:
    topology_prior = build_uniform_rooted_tree_topology_prior(["A", "B", "C", "D"])
    branch_length_prior = build_fixed_branch_length_prior(fixed_value=0.5)
    substitution_prior = build_substitution_parameter_prior_bundle(
        kappa_prior=build_fixed_positive_substitution_parameter_prior(fixed_value=2.0),
        exchangeability_prior=build_fixed_simplex_substitution_parameter_prior(
            expected_component_names=DNA_EXCHANGEABILITY_LABELS,
            fixed_values={
                "AC": 0.1,
                "AG": 0.2,
                "AT": 0.15,
                "CG": 0.1,
                "CT": 0.25,
                "GT": 0.2,
            },
        ),
        base_frequency_prior=build_fixed_simplex_substitution_parameter_prior(
            expected_component_names=("A", "C", "G", "T"),
            fixed_values={
                "A": 0.3,
                "C": 0.2,
                "G": 0.1,
                "T": 0.4,
            },
        ),
        gamma_alpha_prior=build_fixed_positive_substitution_parameter_prior(
            fixed_value=1.5
        ),
    )

    sample = sample_prior_only_phylogenetic_state(
        tree_topology_prior=topology_prior,
        branch_length_prior=branch_length_prior,
        substitution_parameter_prior=substitution_prior,
        seed=5,
    )

    assert sample.substitution_parameter_state == PriorOnlySubstitutionParameterState(
        kappa=2.0,
        exchangeabilities={
            "AC": 0.1,
            "AG": 0.2,
            "AT": 0.15,
            "CG": 0.1,
            "CT": 0.25,
            "GT": 0.2,
        },
        base_frequencies={
            "A": 0.3,
            "C": 0.2,
            "G": 0.1,
            "T": 0.4,
        },
        gamma_alpha=1.5,
        invariant_proportion=None,
    )
    assert all(
        math.isclose(
            branch_row.branch_length,
            0.5,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for branch_row in sample.branch_rows
    )
    assert math.isclose(sample.branch_length_log_prior, 0.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(sample.substitution_log_prior, 0.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        sample.total_log_prior,
        sample.topology_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_prior_only_simulation_requires_positive_sample_count() -> None:
    topology_prior = build_uniform_rooted_tree_topology_prior(["A", "B"])
    branch_length_prior = build_exponential_branch_length_prior(rate=1.0)

    with pytest.raises(
        PhylogeneticsError,
        match="requires sample_count to be positive",
    ):
        simulate_prior_only_phylogenetic_states(
            tree_topology_prior=topology_prior,
            branch_length_prior=branch_length_prior,
            sample_count=0,
            seed=0,
        )
