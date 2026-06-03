from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    FIXED_TOPOLOGY_STRICT_CLOCK_MODELS,
    FixedTopologyStrictClockModelDefinition,
    FixedTopologyStrictClockNodeAgeSummary,
    FixedTopologyStrictClockPosteriorRow,
    FixedTopologyStrictClockProposalSchedule,
    FixedTopologyStrictClockRateSummary,
    FixedTopologyStrictClockRunReport,
    build_fixed_topology_strict_clock_model_definition,
    build_fixed_topology_strict_clock_proposal_schedule,
    run_fixed_topology_strict_clock_metropolis_hastings,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    FIXED_TOPOLOGY_STRICT_CLOCK_MODELS as FIXED_TOPOLOGY_STRICT_CLOCK_MODELS_IMPL,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    FixedTopologyStrictClockModelDefinition as FixedTopologyStrictClockModelDefinitionImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    FixedTopologyStrictClockNodeAgeSummary as FixedTopologyStrictClockNodeAgeSummaryImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    FixedTopologyStrictClockPosteriorRow as FixedTopologyStrictClockPosteriorRowImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    FixedTopologyStrictClockProposalSchedule as FixedTopologyStrictClockProposalScheduleImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    FixedTopologyStrictClockRateSummary as FixedTopologyStrictClockRateSummaryImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    FixedTopologyStrictClockRunReport as FixedTopologyStrictClockRunReportImpl,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    build_fixed_topology_strict_clock_model_definition as build_fixed_topology_strict_clock_model_definition_impl,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    build_fixed_topology_strict_clock_proposal_schedule as build_fixed_topology_strict_clock_proposal_schedule_impl,
)
from bijux_phylogenetics.bayesian.fixed_topology_strict_clock import (
    run_fixed_topology_strict_clock_metropolis_hastings as run_fixed_topology_strict_clock_metropolis_hastings_impl,
)


def test_bayesian_exports_fixed_topology_strict_clock_surface() -> None:
    assert FIXED_TOPOLOGY_STRICT_CLOCK_MODELS == FIXED_TOPOLOGY_STRICT_CLOCK_MODELS_IMPL
    assert (
        FixedTopologyStrictClockModelDefinition
        is FixedTopologyStrictClockModelDefinitionImpl
    )
    assert (
        FixedTopologyStrictClockNodeAgeSummary
        is FixedTopologyStrictClockNodeAgeSummaryImpl
    )
    assert (
        FixedTopologyStrictClockPosteriorRow is FixedTopologyStrictClockPosteriorRowImpl
    )
    assert (
        FixedTopologyStrictClockProposalSchedule
        is FixedTopologyStrictClockProposalScheduleImpl
    )
    assert (
        FixedTopologyStrictClockRateSummary is FixedTopologyStrictClockRateSummaryImpl
    )
    assert FixedTopologyStrictClockRunReport is FixedTopologyStrictClockRunReportImpl
    assert (
        build_fixed_topology_strict_clock_model_definition
        is build_fixed_topology_strict_clock_model_definition_impl
    )
    assert (
        build_fixed_topology_strict_clock_proposal_schedule
        is build_fixed_topology_strict_clock_proposal_schedule_impl
    )
    assert (
        run_fixed_topology_strict_clock_metropolis_hastings
        is run_fixed_topology_strict_clock_metropolis_hastings_impl
    )
