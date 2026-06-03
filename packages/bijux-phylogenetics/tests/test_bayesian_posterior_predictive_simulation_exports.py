from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES,
    PosteriorPredictiveAlignmentReplicate,
    PosteriorPredictiveAlignmentSimulationReport,
    PosteriorPredictiveContinuousTraitReplicate,
    PosteriorPredictiveContinuousTraitSimulationReport,
    PosteriorPredictiveDiscreteTraitReplicate,
    PosteriorPredictiveDiscreteTraitSimulationReport,
    PosteriorPredictiveObservedStatisticRow,
    PosteriorPredictiveReplicateStatisticRow,
    PosteriorPredictiveSimulationDefinition,
    PosteriorPredictiveStatisticSummaryRow,
    build_posterior_predictive_simulation_definition,
    simulate_brownian_continuous_trait_posterior_predictive,
    simulate_discrete_trait_mk_posterior_predictive,
    simulate_fixed_topology_dna_posterior_predictive,
    simulate_fixed_topology_partitioned_dna_posterior_predictive,
    simulate_joint_topology_dna_posterior_predictive,
    simulate_ornstein_uhlenbeck_continuous_trait_posterior_predictive,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES as POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES_IMPL,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    PosteriorPredictiveAlignmentReplicate as PosteriorPredictiveAlignmentReplicateImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    PosteriorPredictiveAlignmentSimulationReport as PosteriorPredictiveAlignmentSimulationReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    PosteriorPredictiveContinuousTraitReplicate as PosteriorPredictiveContinuousTraitReplicateImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    PosteriorPredictiveContinuousTraitSimulationReport as PosteriorPredictiveContinuousTraitSimulationReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    PosteriorPredictiveDiscreteTraitReplicate as PosteriorPredictiveDiscreteTraitReplicateImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    PosteriorPredictiveDiscreteTraitSimulationReport as PosteriorPredictiveDiscreteTraitSimulationReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    PosteriorPredictiveObservedStatisticRow as PosteriorPredictiveObservedStatisticRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    PosteriorPredictiveReplicateStatisticRow as PosteriorPredictiveReplicateStatisticRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    PosteriorPredictiveSimulationDefinition as PosteriorPredictiveSimulationDefinitionImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    PosteriorPredictiveStatisticSummaryRow as PosteriorPredictiveStatisticSummaryRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    build_posterior_predictive_simulation_definition as build_posterior_predictive_simulation_definition_impl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    simulate_brownian_continuous_trait_posterior_predictive as simulate_brownian_continuous_trait_posterior_predictive_impl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    simulate_discrete_trait_mk_posterior_predictive as simulate_discrete_trait_mk_posterior_predictive_impl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    simulate_fixed_topology_dna_posterior_predictive as simulate_fixed_topology_dna_posterior_predictive_impl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    simulate_fixed_topology_partitioned_dna_posterior_predictive as simulate_fixed_topology_partitioned_dna_posterior_predictive_impl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    simulate_joint_topology_dna_posterior_predictive as simulate_joint_topology_dna_posterior_predictive_impl,
)
from bijux_phylogenetics.bayesian.posterior_predictive_simulation import (
    simulate_ornstein_uhlenbeck_continuous_trait_posterior_predictive as simulate_ornstein_uhlenbeck_continuous_trait_posterior_predictive_impl,
)


def test_bayesian_exports_posterior_predictive_simulation_surface() -> None:
    assert (
        POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES
        is POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES_IMPL
    )
    assert (
        PosteriorPredictiveSimulationDefinition
        is PosteriorPredictiveSimulationDefinitionImpl
    )
    assert (
        PosteriorPredictiveObservedStatisticRow
        is PosteriorPredictiveObservedStatisticRowImpl
    )
    assert (
        PosteriorPredictiveReplicateStatisticRow
        is PosteriorPredictiveReplicateStatisticRowImpl
    )
    assert (
        PosteriorPredictiveStatisticSummaryRow
        is PosteriorPredictiveStatisticSummaryRowImpl
    )
    assert (
        PosteriorPredictiveAlignmentReplicate
        is PosteriorPredictiveAlignmentReplicateImpl
    )
    assert (
        PosteriorPredictiveDiscreteTraitReplicate
        is PosteriorPredictiveDiscreteTraitReplicateImpl
    )
    assert (
        PosteriorPredictiveContinuousTraitReplicate
        is PosteriorPredictiveContinuousTraitReplicateImpl
    )
    assert (
        PosteriorPredictiveAlignmentSimulationReport
        is PosteriorPredictiveAlignmentSimulationReportImpl
    )
    assert (
        PosteriorPredictiveDiscreteTraitSimulationReport
        is PosteriorPredictiveDiscreteTraitSimulationReportImpl
    )
    assert (
        PosteriorPredictiveContinuousTraitSimulationReport
        is PosteriorPredictiveContinuousTraitSimulationReportImpl
    )
    assert (
        build_posterior_predictive_simulation_definition
        is build_posterior_predictive_simulation_definition_impl
    )
    assert (
        simulate_fixed_topology_dna_posterior_predictive
        is simulate_fixed_topology_dna_posterior_predictive_impl
    )
    assert (
        simulate_joint_topology_dna_posterior_predictive
        is simulate_joint_topology_dna_posterior_predictive_impl
    )
    assert (
        simulate_fixed_topology_partitioned_dna_posterior_predictive
        is simulate_fixed_topology_partitioned_dna_posterior_predictive_impl
    )
    assert (
        simulate_discrete_trait_mk_posterior_predictive
        is simulate_discrete_trait_mk_posterior_predictive_impl
    )
    assert (
        simulate_brownian_continuous_trait_posterior_predictive
        is simulate_brownian_continuous_trait_posterior_predictive_impl
    )
    assert (
        simulate_ornstein_uhlenbeck_continuous_trait_posterior_predictive
        is simulate_ornstein_uhlenbeck_continuous_trait_posterior_predictive_impl
    )
