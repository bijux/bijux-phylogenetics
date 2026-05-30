from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaModelDefinition,
    FixedTopologyDnaProposalSchedule,
)
from bijux_phylogenetics.bayesian.tree_topology_priors import TreeTopologyPriorModel
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

JOINT_TOPOLOGY_DNA_TOPOLOGY_MOVE_KINDS = ("nni", "spr", "tbr")


@dataclass(frozen=True, slots=True)
class JointTopologyDnaModelDefinition:
    """One validated joint topology-and-parameter DNA posterior model."""

    sequence_model_definition: FixedTopologyDnaModelDefinition
    topology_prior: TreeTopologyPriorModel

    @property
    def substitution_model_name(self) -> str:
        return self.sequence_model_definition.substitution_model_name

    @property
    def active_parameter_targets(self) -> tuple[str, ...]:
        return self.sequence_model_definition.active_parameter_targets


@dataclass(frozen=True, slots=True)
class JointTopologyDnaProposalSchedule:
    """One validated proposal schedule for joint DNA topology and parameter sampling."""

    sequence_proposal_schedule: FixedTopologyDnaProposalSchedule
    nni_move_weight: float
    spr_move_weight: float
    tbr_move_weight: float

    @property
    def substitution_model_name(self) -> str:
        return self.sequence_proposal_schedule.substitution_model_name


def build_joint_topology_dna_model_definition(
    *,
    sequence_model_definition: FixedTopologyDnaModelDefinition,
    topology_prior: TreeTopologyPriorModel,
) -> JointTopologyDnaModelDefinition:
    """Build one validated joint DNA topology-and-parameter posterior model."""
    if not isinstance(sequence_model_definition, FixedTopologyDnaModelDefinition):
        raise PhylogeneticsError(
            "joint topology DNA posterior model requires one FixedTopologyDnaModelDefinition",
            code="joint_topology_dna_sequence_model_definition_type_invalid",
        )
    if not isinstance(topology_prior, TreeTopologyPriorModel):
        raise PhylogeneticsError(
            "joint topology DNA posterior model requires one TreeTopologyPriorModel",
            code="joint_topology_dna_topology_prior_type_invalid",
        )
    return JointTopologyDnaModelDefinition(
        sequence_model_definition=sequence_model_definition,
        topology_prior=topology_prior,
    )


def build_joint_topology_dna_proposal_schedule(
    *,
    sequence_proposal_schedule: FixedTopologyDnaProposalSchedule,
    nni_move_weight: float = 0.0,
    spr_move_weight: float = 0.0,
    tbr_move_weight: float = 0.0,
) -> JointTopologyDnaProposalSchedule:
    """Build one validated joint DNA proposal schedule with explicit topology moves."""
    if not isinstance(sequence_proposal_schedule, FixedTopologyDnaProposalSchedule):
        raise PhylogeneticsError(
            "joint topology DNA proposal schedule requires one FixedTopologyDnaProposalSchedule",
            code="joint_topology_dna_sequence_proposal_schedule_type_invalid",
        )
    validated_nni_move_weight = _validate_nonnegative_finite_float(
        value=nni_move_weight,
        field_name="nni_move_weight",
        owner_name="joint topology DNA proposal schedule",
    )
    validated_spr_move_weight = _validate_nonnegative_finite_float(
        value=spr_move_weight,
        field_name="spr_move_weight",
        owner_name="joint topology DNA proposal schedule",
    )
    validated_tbr_move_weight = _validate_nonnegative_finite_float(
        value=tbr_move_weight,
        field_name="tbr_move_weight",
        owner_name="joint topology DNA proposal schedule",
    )
    topology_weight_total = math.fsum(
        (
            validated_nni_move_weight,
            validated_spr_move_weight,
            validated_tbr_move_weight,
        )
    )
    if topology_weight_total <= 0.0:
        raise PhylogeneticsError(
            "joint topology DNA proposal schedule requires at least one positive topology move weight",
            code="joint_topology_dna_topology_move_weight_missing",
        )
    return JointTopologyDnaProposalSchedule(
        sequence_proposal_schedule=sequence_proposal_schedule,
        nni_move_weight=validated_nni_move_weight,
        spr_move_weight=validated_spr_move_weight,
        tbr_move_weight=validated_tbr_move_weight,
    )


def _validate_nonnegative_finite_float(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    try:
        validated_value = float(value)
    except (TypeError, ValueError) as error:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be a finite float",
            code="joint_topology_dna_float_type_invalid",
            details={"field_name": field_name},
        ) from error
    if not math.isfinite(validated_value) or validated_value < 0.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be nonnegative and finite",
            code="joint_topology_dna_nonnegative_float_invalid",
            details={"field_name": field_name, "value": value},
        )
    return float(format(validated_value, ".15g"))


__all__ = [
    "JOINT_TOPOLOGY_DNA_TOPOLOGY_MOVE_KINDS",
    "JointTopologyDnaModelDefinition",
    "JointTopologyDnaProposalSchedule",
    "build_joint_topology_dna_model_definition",
    "build_joint_topology_dna_proposal_schedule",
]
