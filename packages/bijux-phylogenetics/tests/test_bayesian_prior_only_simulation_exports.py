from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    PriorOnlyPhylogeneticSample,
    PriorOnlyPhylogeneticSimulationReport,
    PriorOnlySampledBranchRow,
    PriorOnlySubstitutionParameterState,
    sample_prior_only_phylogenetic_state,
    simulate_prior_only_phylogenetic_states,
)
from bijux_phylogenetics.bayesian.prior_sampling import (
    PriorOnlyPhylogeneticSample as PriorOnlyPhylogeneticSampleImpl,
)
from bijux_phylogenetics.bayesian.prior_sampling import (
    PriorOnlyPhylogeneticSimulationReport as PriorOnlyPhylogeneticSimulationReportImpl,
)
from bijux_phylogenetics.bayesian.prior_sampling import (
    PriorOnlySampledBranchRow as PriorOnlySampledBranchRowImpl,
)
from bijux_phylogenetics.bayesian.prior_sampling import (
    PriorOnlySubstitutionParameterState as PriorOnlySubstitutionParameterStateImpl,
)
from bijux_phylogenetics.bayesian.prior_sampling import (
    sample_prior_only_phylogenetic_state as sample_prior_only_phylogenetic_state_impl,
)
from bijux_phylogenetics.bayesian.prior_sampling import (
    simulate_prior_only_phylogenetic_states as simulate_prior_only_phylogenetic_states_impl,
)


def test_bayesian_exports_prior_only_simulation_surface() -> None:
    assert PriorOnlySampledBranchRow is PriorOnlySampledBranchRowImpl
    assert (
        PriorOnlySubstitutionParameterState is PriorOnlySubstitutionParameterStateImpl
    )
    assert PriorOnlyPhylogeneticSample is PriorOnlyPhylogeneticSampleImpl
    assert (
        PriorOnlyPhylogeneticSimulationReport
        is PriorOnlyPhylogeneticSimulationReportImpl
    )
    assert (
        sample_prior_only_phylogenetic_state
        is sample_prior_only_phylogenetic_state_impl
    )
    assert (
        simulate_prior_only_phylogenetic_states
        is simulate_prior_only_phylogenetic_states_impl
    )
