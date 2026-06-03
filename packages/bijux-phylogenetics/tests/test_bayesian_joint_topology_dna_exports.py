from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    JOINT_TOPOLOGY_DNA_TOPOLOGY_MOVE_KINDS,
    JointTopologyDnaModelDefinition,
    JointTopologyDnaPosteriorRow,
    JointTopologyDnaProposalSchedule,
    JointTopologyDnaRunReport,
    build_joint_topology_dna_model_definition,
    build_joint_topology_dna_proposal_schedule,
    run_joint_topology_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    JOINT_TOPOLOGY_DNA_TOPOLOGY_MOVE_KINDS as JOINT_TOPOLOGY_DNA_TOPOLOGY_MOVE_KINDS_IMPL,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    JointTopologyDnaModelDefinition as JointTopologyDnaModelDefinitionImpl,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    JointTopologyDnaPosteriorRow as JointTopologyDnaPosteriorRowImpl,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    JointTopologyDnaProposalSchedule as JointTopologyDnaProposalScheduleImpl,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    JointTopologyDnaRunReport as JointTopologyDnaRunReportImpl,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    build_joint_topology_dna_model_definition as build_joint_topology_dna_model_definition_impl,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    build_joint_topology_dna_proposal_schedule as build_joint_topology_dna_proposal_schedule_impl,
)
from bijux_phylogenetics.bayesian.joint_topology_dna import (
    run_joint_topology_dna_metropolis_hastings as run_joint_topology_dna_metropolis_hastings_impl,
)


def test_bayesian_exports_joint_topology_dna_surface() -> None:
    assert (
        JOINT_TOPOLOGY_DNA_TOPOLOGY_MOVE_KINDS
        is JOINT_TOPOLOGY_DNA_TOPOLOGY_MOVE_KINDS_IMPL
    )
    assert JointTopologyDnaModelDefinition is JointTopologyDnaModelDefinitionImpl
    assert JointTopologyDnaPosteriorRow is JointTopologyDnaPosteriorRowImpl
    assert JointTopologyDnaProposalSchedule is JointTopologyDnaProposalScheduleImpl
    assert JointTopologyDnaRunReport is JointTopologyDnaRunReportImpl
    assert (
        build_joint_topology_dna_model_definition
        is build_joint_topology_dna_model_definition_impl
    )
    assert (
        build_joint_topology_dna_proposal_schedule
        is build_joint_topology_dna_proposal_schedule_impl
    )
    assert (
        run_joint_topology_dna_metropolis_hastings
        is run_joint_topology_dna_metropolis_hastings_impl
    )
