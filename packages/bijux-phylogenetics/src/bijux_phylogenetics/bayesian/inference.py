from __future__ import annotations

from collections.abc import Sequence

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .fixed_topology_dna import (
    FixedTopologyDnaModelDefinition,
    FixedTopologyDnaProposalSchedule,
    FixedTopologyDnaRunReport,
    run_fixed_topology_dna_metropolis_hastings,
)

__all__ = ["run_bayesian_inference"]


def run_bayesian_inference(
    *,
    tree: PhyloTree,
    records: Sequence[AlignmentRecord],
    model_definition: FixedTopologyDnaModelDefinition,
    proposal_schedule: FixedTopologyDnaProposalSchedule,
    iteration_count: int,
    sample_every: int = 1,
    seed: int = 0,
    observation_policy: str = "reject",
) -> FixedTopologyDnaRunReport:
    """Run one native Bayesian inference chain through the stable public entry point."""
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
