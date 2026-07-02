from __future__ import annotations

from collections.abc import Sequence
from typing import TypeAlias

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .fixed_topology_dna import (
    FixedTopologyDnaModelDefinition,
    FixedTopologyDnaProposalSchedule,
    FixedTopologyDnaRunReport,
    run_fixed_topology_dna_metropolis_hastings,
)
from .joint_topology_dna import (
    JointTopologyDnaModelDefinition,
    JointTopologyDnaProposalSchedule,
    JointTopologyDnaRunReport,
    run_joint_topology_dna_metropolis_hastings,
)

BayesianInferenceModelDefinition: TypeAlias = (
    FixedTopologyDnaModelDefinition | JointTopologyDnaModelDefinition
)
BayesianInferenceProposalSchedule: TypeAlias = (
    FixedTopologyDnaProposalSchedule | JointTopologyDnaProposalSchedule
)
BayesianInferenceRunReport: TypeAlias = (
    FixedTopologyDnaRunReport | JointTopologyDnaRunReport
)

__all__ = [
    "BayesianInferenceModelDefinition",
    "BayesianInferenceProposalSchedule",
    "BayesianInferenceRunReport",
    "run_bayesian_inference",
]


def run_bayesian_inference(
    *,
    tree: PhyloTree,
    records: Sequence[AlignmentRecord],
    model_definition: BayesianInferenceModelDefinition,
    proposal_schedule: BayesianInferenceProposalSchedule,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
    observation_policy: str = "reject",
) -> BayesianInferenceRunReport:
    """Run one native Bayesian inference chain through the stable public entry point."""
    if isinstance(model_definition, FixedTopologyDnaModelDefinition) and isinstance(
        proposal_schedule, FixedTopologyDnaProposalSchedule
    ):
        return run_fixed_topology_dna_metropolis_hastings(
            tree=tree,
            records=records,
            model_definition=model_definition,
            proposal_schedule=proposal_schedule,
            iteration_count=iteration_count,
            sample_every=sample_every,
            seed=seed,
            observation_policy=observation_policy,
        )
    if isinstance(model_definition, JointTopologyDnaModelDefinition) and isinstance(
        proposal_schedule, JointTopologyDnaProposalSchedule
    ):
        return run_joint_topology_dna_metropolis_hastings(
            tree=tree,
            records=records,
            model_definition=model_definition,
            proposal_schedule=proposal_schedule,
            iteration_count=iteration_count,
            sample_every=sample_every,
            seed=seed,
            observation_policy=observation_policy,
        )
    raise PhylogeneticsError(
        "native Bayesian inference requires one supported model-definition and proposal-schedule pair",
        code="bayesian_inference_model_schedule_unsupported",
        details={
            "model_definition_type": type(model_definition).__name__,
            "proposal_schedule_type": type(proposal_schedule).__name__,
        },
    )
