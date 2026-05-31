from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
    build_fixed_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaModelDefinition,
    FixedTopologyDnaProposalSchedule,
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_dirichlet_simplex_substitution_parameter_prior,
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_build_fixed_topology_dna_model_definition_records_active_parameter_surface() -> (
    None
):
    model_definition = build_fixed_topology_dna_model_definition(
        substitution_model_name="HKY85",
        branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
        substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=1.5
            ),
            base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                expected_component_names=("A", "C", "G", "T"),
                concentration_parameters={"A": 2.0, "C": 2.0, "G": 2.0, "T": 2.0},
            ),
        ),
        initial_kappa=2.5,
        initial_base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
    )

    assert isinstance(model_definition, FixedTopologyDnaModelDefinition)
    assert model_definition.substitution_model_name == "HKY85"
    assert model_definition.active_parameter_targets == ("base-frequencies", "kappa")
    assert model_definition.initial_kappa == 2.5
    assert model_definition.initial_base_frequencies == {
        "A": 0.4,
        "C": 0.1,
        "G": 0.2,
        "T": 0.3,
    }


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        (
            {
                "substitution_model_name": "JC69",
                "branch_length_prior": build_exponential_branch_length_prior(rate=4.0),
                "substitution_parameter_prior_bundle": build_substitution_parameter_prior_bundle(
                    kappa_prior=build_exponential_positive_substitution_parameter_prior(
                        rate=1.5
                    )
                ),
            },
            "does not support JC69",
        ),
        (
            {
                "substitution_model_name": "K80",
                "branch_length_prior": build_fixed_branch_length_prior(fixed_value=0.1),
                "substitution_parameter_prior_bundle": build_substitution_parameter_prior_bundle(
                    kappa_prior=build_exponential_positive_substitution_parameter_prior(
                        rate=1.5
                    )
                ),
            },
            "does not support fixed branch-length priors",
        ),
        (
            {
                "substitution_model_name": "K80",
                "branch_length_prior": build_exponential_branch_length_prior(rate=4.0),
                "substitution_parameter_prior_bundle": build_substitution_parameter_prior_bundle(
                    base_frequency_prior=build_dirichlet_simplex_substitution_parameter_prior(
                        expected_component_names=("A", "C", "G", "T"),
                        concentration_parameters={
                            "A": 2.0,
                            "C": 2.0,
                            "G": 2.0,
                            "T": 2.0,
                        },
                    )
                ),
            },
            "does not sample",
        ),
    ],
)
def test_build_fixed_topology_dna_model_definition_rejects_incompatible_surfaces(
    kwargs: dict[str, object],
    match: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=match):
        build_fixed_topology_dna_model_definition(**kwargs)


def test_build_fixed_topology_dna_proposal_schedule_requires_active_parameter_moves() -> (
    None
):
    model_definition = build_fixed_topology_dna_model_definition(
        substitution_model_name="K80",
        branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
        substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=1.5
            )
        ),
    )

    schedule = build_fixed_topology_dna_proposal_schedule(
        model_definition=model_definition,
        branch_length_move_weight=2.0,
        branch_length_log_scale_standard_deviation=0.3,
        kappa_move_weight=1.0,
        kappa_log_scale_standard_deviation=0.4,
    )

    assert isinstance(schedule, FixedTopologyDnaProposalSchedule)
    assert schedule.substitution_model_name == "K80"
    assert schedule.branch_length_move_weight == 2.0
    assert schedule.kappa_move_weight == 1.0

    with pytest.raises(PhylogeneticsError, match="positive move weight"):
        build_fixed_topology_dna_proposal_schedule(
            model_definition=model_definition,
            branch_length_move_weight=1.0,
            branch_length_log_scale_standard_deviation=0.3,
            kappa_move_weight=0.0,
        )
