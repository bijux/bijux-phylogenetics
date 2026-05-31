from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaModelDefinition,
    FixedTopologyDnaProposalSchedule,
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    JointTopologyDnaModelDefinition,
    JointTopologyDnaProposalSchedule,
    build_joint_topology_dna_model_definition,
    build_joint_topology_dna_proposal_schedule,
)
from bijux_phylogenetics.bayesian.substitution_parameter_priors import (
    build_exponential_positive_substitution_parameter_prior,
    build_substitution_parameter_prior_bundle,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import (
    build_uniform_rooted_tree_topology_prior,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def _build_k80_sequence_model_definition() -> FixedTopologyDnaModelDefinition:
    return build_fixed_topology_dna_model_definition(
        substitution_model_name="K80",
        branch_length_prior=build_exponential_branch_length_prior(rate=4.0),
        substitution_parameter_prior_bundle=build_substitution_parameter_prior_bundle(
            kappa_prior=build_exponential_positive_substitution_parameter_prior(
                rate=1.5
            )
        ),
    )


def _build_k80_sequence_proposal_schedule() -> FixedTopologyDnaProposalSchedule:
    return build_fixed_topology_dna_proposal_schedule(
        model_definition=_build_k80_sequence_model_definition(),
        branch_length_move_weight=1.0,
        branch_length_log_scale_standard_deviation=0.25,
        kappa_move_weight=1.0,
        kappa_log_scale_standard_deviation=0.35,
    )


def test_build_joint_topology_dna_model_definition_composes_sequence_and_topology_surfaces() -> (
    None
):
    sequence_model_definition = _build_k80_sequence_model_definition()
    topology_prior = build_uniform_rooted_tree_topology_prior(["A", "B", "C", "D"])

    model_definition = build_joint_topology_dna_model_definition(
        sequence_model_definition=sequence_model_definition,
        topology_prior=topology_prior,
    )

    assert isinstance(model_definition, JointTopologyDnaModelDefinition)
    assert model_definition.sequence_model_definition is sequence_model_definition
    assert model_definition.topology_prior is topology_prior
    assert model_definition.substitution_model_name == "K80"
    assert model_definition.active_parameter_targets == ("kappa",)


def test_build_joint_topology_dna_proposal_schedule_records_topology_policy() -> None:
    sequence_proposal_schedule = _build_k80_sequence_proposal_schedule()

    proposal_schedule = build_joint_topology_dna_proposal_schedule(
        sequence_proposal_schedule=sequence_proposal_schedule,
        nni_move_weight=2.0,
        spr_move_weight=1.0,
    )

    assert isinstance(proposal_schedule, JointTopologyDnaProposalSchedule)
    assert proposal_schedule.sequence_proposal_schedule is sequence_proposal_schedule
    assert proposal_schedule.substitution_model_name == "K80"
    assert proposal_schedule.nni_move_weight == 2.0
    assert proposal_schedule.spr_move_weight == 1.0
    assert proposal_schedule.tbr_move_weight == 0.0


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        (
            {
                "sequence_model_definition": object(),
                "topology_prior": build_uniform_rooted_tree_topology_prior(
                    ["A", "B", "C", "D"]
                ),
            },
            "requires one FixedTopologyDnaModelDefinition",
        ),
        (
            {
                "sequence_model_definition": _build_k80_sequence_model_definition(),
                "topology_prior": object(),
            },
            "requires one TreeTopologyPriorModel",
        ),
    ],
)
def test_build_joint_topology_dna_model_definition_rejects_invalid_inputs(
    kwargs: dict[str, object],
    match: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=match):
        build_joint_topology_dna_model_definition(**kwargs)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        (
            {
                "sequence_proposal_schedule": object(),
                "nni_move_weight": 1.0,
            },
            "requires one FixedTopologyDnaProposalSchedule",
        ),
        (
            {
                "sequence_proposal_schedule": _build_k80_sequence_proposal_schedule(),
                "nni_move_weight": 0.0,
                "spr_move_weight": 0.0,
                "tbr_move_weight": 0.0,
            },
            "requires at least one positive topology move weight",
        ),
    ],
)
def test_build_joint_topology_dna_proposal_schedule_rejects_invalid_inputs(
    kwargs: dict[str, object],
    match: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=match):
        build_joint_topology_dna_proposal_schedule(**kwargs)
