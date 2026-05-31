from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    FIXED_TOPOLOGY_RELAXED_CLOCK_MODELS,
    FixedTopologyRelaxedClockBranchRateSummary,
    FixedTopologyRelaxedClockModelDefinition,
    FixedTopologyRelaxedClockNodeAgeSummary,
    FixedTopologyRelaxedClockPosteriorRow,
    FixedTopologyRelaxedClockProposalSchedule,
    FixedTopologyRelaxedClockRunReport,
    build_fixed_topology_relaxed_clock_model_definition,
    build_fixed_topology_relaxed_clock_proposal_schedule,
    run_fixed_topology_relaxed_clock_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    FIXED_TOPOLOGY_RELAXED_CLOCK_MODELS as FIXED_TOPOLOGY_RELAXED_CLOCK_MODELS_IMPL,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    FixedTopologyRelaxedClockBranchRateSummary as FixedTopologyRelaxedClockBranchRateSummaryImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    FixedTopologyRelaxedClockModelDefinition as FixedTopologyRelaxedClockModelDefinitionImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    FixedTopologyRelaxedClockNodeAgeSummary as FixedTopologyRelaxedClockNodeAgeSummaryImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    FixedTopologyRelaxedClockPosteriorRow as FixedTopologyRelaxedClockPosteriorRowImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    FixedTopologyRelaxedClockProposalSchedule as FixedTopologyRelaxedClockProposalScheduleImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    FixedTopologyRelaxedClockRunReport as FixedTopologyRelaxedClockRunReportImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    build_fixed_topology_relaxed_clock_model_definition as build_fixed_topology_relaxed_clock_model_definition_impl,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    build_fixed_topology_relaxed_clock_proposal_schedule as build_fixed_topology_relaxed_clock_proposal_schedule_impl,
)
from bijux_phylogenetics.bayesian.fixed_topology_relaxed_clock import (
    run_fixed_topology_relaxed_clock_metropolis_hastings as run_fixed_topology_relaxed_clock_metropolis_hastings_impl,
)


def test_bayesian_exports_fixed_topology_relaxed_clock_surface() -> None:
    assert (
        FIXED_TOPOLOGY_RELAXED_CLOCK_MODELS == FIXED_TOPOLOGY_RELAXED_CLOCK_MODELS_IMPL
    )
    assert (
        FixedTopologyRelaxedClockBranchRateSummary
        is FixedTopologyRelaxedClockBranchRateSummaryImpl
    )
    assert (
        FixedTopologyRelaxedClockModelDefinition
        is FixedTopologyRelaxedClockModelDefinitionImpl
    )
    assert (
        FixedTopologyRelaxedClockNodeAgeSummary
        is FixedTopologyRelaxedClockNodeAgeSummaryImpl
    )
    assert (
        FixedTopologyRelaxedClockPosteriorRow
        is FixedTopologyRelaxedClockPosteriorRowImpl
    )
    assert (
        FixedTopologyRelaxedClockProposalSchedule
        is FixedTopologyRelaxedClockProposalScheduleImpl
    )
    assert FixedTopologyRelaxedClockRunReport is FixedTopologyRelaxedClockRunReportImpl
    assert (
        build_fixed_topology_relaxed_clock_model_definition
        is build_fixed_topology_relaxed_clock_model_definition_impl
    )
    assert (
        build_fixed_topology_relaxed_clock_proposal_schedule
        is build_fixed_topology_relaxed_clock_proposal_schedule_impl
    )
    assert (
        run_fixed_topology_relaxed_clock_metropolis_hastings
        is run_fixed_topology_relaxed_clock_metropolis_hastings_impl
    )
