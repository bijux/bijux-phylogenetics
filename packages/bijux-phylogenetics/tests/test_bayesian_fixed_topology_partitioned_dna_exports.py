from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    FIXED_TOPOLOGY_PARTITIONED_DNA_SUBSTITUTION_MODELS,
    FixedTopologyPartitionedDnaModelDefinition,
    FixedTopologyPartitionedDnaPartitionRow,
    FixedTopologyPartitionedDnaPosteriorRow,
    FixedTopologyPartitionedDnaProposalSchedule,
    FixedTopologyPartitionedDnaRunReport,
    build_fixed_topology_partitioned_dna_model_definition,
    build_fixed_topology_partitioned_dna_proposal_schedule,
    run_fixed_topology_partitioned_dna_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    FIXED_TOPOLOGY_PARTITIONED_DNA_SUBSTITUTION_MODELS as FIXED_TOPOLOGY_PARTITIONED_DNA_SUBSTITUTION_MODELS_IMPL,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    FixedTopologyPartitionedDnaModelDefinition as FixedTopologyPartitionedDnaModelDefinitionImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    FixedTopologyPartitionedDnaPartitionRow as FixedTopologyPartitionedDnaPartitionRowImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    FixedTopologyPartitionedDnaPosteriorRow as FixedTopologyPartitionedDnaPosteriorRowImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    FixedTopologyPartitionedDnaProposalSchedule as FixedTopologyPartitionedDnaProposalScheduleImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    FixedTopologyPartitionedDnaRunReport as FixedTopologyPartitionedDnaRunReportImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    build_fixed_topology_partitioned_dna_model_definition as build_fixed_topology_partitioned_dna_model_definition_impl,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    build_fixed_topology_partitioned_dna_proposal_schedule as build_fixed_topology_partitioned_dna_proposal_schedule_impl,
)
from bijux_phylogenetics.bayesian.fixed_topology_partitioned_dna import (
    run_fixed_topology_partitioned_dna_metropolis_hastings as run_fixed_topology_partitioned_dna_metropolis_hastings_impl,
)


def test_bayesian_exports_fixed_topology_partitioned_dna_surface() -> None:
    assert (
        FIXED_TOPOLOGY_PARTITIONED_DNA_SUBSTITUTION_MODELS
        is FIXED_TOPOLOGY_PARTITIONED_DNA_SUBSTITUTION_MODELS_IMPL
    )
    assert (
        FixedTopologyPartitionedDnaModelDefinition
        is FixedTopologyPartitionedDnaModelDefinitionImpl
    )
    assert (
        FixedTopologyPartitionedDnaPartitionRow
        is FixedTopologyPartitionedDnaPartitionRowImpl
    )
    assert (
        FixedTopologyPartitionedDnaPosteriorRow
        is FixedTopologyPartitionedDnaPosteriorRowImpl
    )
    assert (
        FixedTopologyPartitionedDnaProposalSchedule
        is FixedTopologyPartitionedDnaProposalScheduleImpl
    )
    assert (
        FixedTopologyPartitionedDnaRunReport is FixedTopologyPartitionedDnaRunReportImpl
    )
    assert (
        build_fixed_topology_partitioned_dna_model_definition
        is build_fixed_topology_partitioned_dna_model_definition_impl
    )
    assert (
        build_fixed_topology_partitioned_dna_proposal_schedule
        is build_fixed_topology_partitioned_dna_proposal_schedule_impl
    )
    assert (
        run_fixed_topology_partitioned_dna_metropolis_hastings
        is run_fixed_topology_partitioned_dna_metropolis_hastings_impl
    )
