from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    FIXED_TOPOLOGY_DNA_SUBSTITUTION_MODELS,
    FixedTopologyDnaModelDefinition,
    FixedTopologyDnaPosteriorRow,
    FixedTopologyDnaProposalSchedule,
    FixedTopologyDnaRunReport,
    build_fixed_topology_dna_model_definition,
    build_fixed_topology_dna_proposal_schedule,
    run_fixed_topology_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FIXED_TOPOLOGY_DNA_SUBSTITUTION_MODELS as FIXED_TOPOLOGY_DNA_SUBSTITUTION_MODELS_IMPL,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaModelDefinition as FixedTopologyDnaModelDefinitionImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaPosteriorRow as FixedTopologyDnaPosteriorRowImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaProposalSchedule as FixedTopologyDnaProposalScheduleImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    FixedTopologyDnaRunReport as FixedTopologyDnaRunReportImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    build_fixed_topology_dna_model_definition as build_fixed_topology_dna_model_definition_impl,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    build_fixed_topology_dna_proposal_schedule as build_fixed_topology_dna_proposal_schedule_impl,
)
from bijux_phylogenetics.bayesian.fixed_topology_dna import (
    run_fixed_topology_dna_metropolis_hastings as run_fixed_topology_dna_metropolis_hastings_impl,
)


def test_bayesian_exports_fixed_topology_dna_surface() -> None:
    assert (
        FIXED_TOPOLOGY_DNA_SUBSTITUTION_MODELS
        is FIXED_TOPOLOGY_DNA_SUBSTITUTION_MODELS_IMPL
    )
    assert FixedTopologyDnaModelDefinition is FixedTopologyDnaModelDefinitionImpl
    assert FixedTopologyDnaPosteriorRow is FixedTopologyDnaPosteriorRowImpl
    assert FixedTopologyDnaProposalSchedule is FixedTopologyDnaProposalScheduleImpl
    assert FixedTopologyDnaRunReport is FixedTopologyDnaRunReportImpl
    assert (
        build_fixed_topology_dna_model_definition
        is build_fixed_topology_dna_model_definition_impl
    )
    assert (
        build_fixed_topology_dna_proposal_schedule
        is build_fixed_topology_dna_proposal_schedule_impl
    )
    assert (
        run_fixed_topology_dna_metropolis_hastings
        is run_fixed_topology_dna_metropolis_hastings_impl
    )
