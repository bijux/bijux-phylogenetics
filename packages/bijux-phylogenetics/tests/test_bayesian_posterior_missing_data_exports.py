from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    PosteriorMissingContinuousTraitDefinition,
    PosteriorMissingContinuousTraitReport,
    PosteriorMissingContinuousTraitTaxonSummaryRow,
    PosteriorMissingDiscreteTraitDefinition,
    PosteriorMissingDiscreteTraitReport,
    PosteriorMissingDiscreteTraitStateProbabilityRow,
    PosteriorMissingDiscreteTraitTaxonSummaryRow,
    PosteriorMissingNucleotideDefinition,
    PosteriorMissingNucleotideReport,
    PosteriorMissingNucleotideSequenceRecord,
    PosteriorMissingNucleotideSiteSummaryRow,
    PosteriorMissingNucleotideStateProbabilityRow,
    build_posterior_missing_continuous_trait_definition,
    build_posterior_missing_discrete_trait_definition,
    build_posterior_missing_nucleotide_definition,
    summarize_brownian_continuous_trait_posterior_missing_values,
    summarize_continuous_trait_posterior_missing_values,
    summarize_discrete_trait_mk_posterior_missing_states,
    summarize_fixed_topology_dna_posterior_missing_states,
    summarize_fixed_topology_partitioned_dna_posterior_missing_states,
    summarize_joint_topology_dna_posterior_missing_states,
    summarize_nucleotide_posterior_missing_states,
    summarize_ornstein_uhlenbeck_continuous_trait_posterior_missing_values,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingContinuousTraitDefinition as PosteriorMissingContinuousTraitDefinitionImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingContinuousTraitReport as PosteriorMissingContinuousTraitReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingContinuousTraitTaxonSummaryRow as PosteriorMissingContinuousTraitTaxonSummaryRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingDiscreteTraitDefinition as PosteriorMissingDiscreteTraitDefinitionImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingDiscreteTraitReport as PosteriorMissingDiscreteTraitReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingDiscreteTraitStateProbabilityRow as PosteriorMissingDiscreteTraitStateProbabilityRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingDiscreteTraitTaxonSummaryRow as PosteriorMissingDiscreteTraitTaxonSummaryRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingNucleotideDefinition as PosteriorMissingNucleotideDefinitionImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingNucleotideReport as PosteriorMissingNucleotideReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingNucleotideSequenceRecord as PosteriorMissingNucleotideSequenceRecordImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingNucleotideSiteSummaryRow as PosteriorMissingNucleotideSiteSummaryRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    PosteriorMissingNucleotideStateProbabilityRow as PosteriorMissingNucleotideStateProbabilityRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    build_posterior_missing_continuous_trait_definition as build_posterior_missing_continuous_trait_definition_impl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    build_posterior_missing_discrete_trait_definition as build_posterior_missing_discrete_trait_definition_impl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    build_posterior_missing_nucleotide_definition as build_posterior_missing_nucleotide_definition_impl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    summarize_brownian_continuous_trait_posterior_missing_values as summarize_brownian_continuous_trait_posterior_missing_values_impl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    summarize_continuous_trait_posterior_missing_values as summarize_continuous_trait_posterior_missing_values_impl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    summarize_discrete_trait_mk_posterior_missing_states as summarize_discrete_trait_mk_posterior_missing_states_impl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    summarize_fixed_topology_dna_posterior_missing_states as summarize_fixed_topology_dna_posterior_missing_states_impl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    summarize_fixed_topology_partitioned_dna_posterior_missing_states as summarize_fixed_topology_partitioned_dna_posterior_missing_states_impl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    summarize_joint_topology_dna_posterior_missing_states as summarize_joint_topology_dna_posterior_missing_states_impl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    summarize_nucleotide_posterior_missing_states as summarize_nucleotide_posterior_missing_states_impl,
)
from bijux_phylogenetics.bayesian.posterior_missing_data import (
    summarize_ornstein_uhlenbeck_continuous_trait_posterior_missing_values as summarize_ornstein_uhlenbeck_continuous_trait_posterior_missing_values_impl,
)


def test_bayesian_exports_posterior_missing_data_surface() -> None:
    assert (
        PosteriorMissingContinuousTraitDefinition
        is PosteriorMissingContinuousTraitDefinitionImpl
    )
    assert (
        PosteriorMissingContinuousTraitReport
        is PosteriorMissingContinuousTraitReportImpl
    )
    assert (
        PosteriorMissingContinuousTraitTaxonSummaryRow
        is PosteriorMissingContinuousTraitTaxonSummaryRowImpl
    )
    assert (
        PosteriorMissingDiscreteTraitDefinition
        is PosteriorMissingDiscreteTraitDefinitionImpl
    )
    assert (
        PosteriorMissingDiscreteTraitReport is PosteriorMissingDiscreteTraitReportImpl
    )
    assert (
        PosteriorMissingDiscreteTraitStateProbabilityRow
        is PosteriorMissingDiscreteTraitStateProbabilityRowImpl
    )
    assert (
        PosteriorMissingDiscreteTraitTaxonSummaryRow
        is PosteriorMissingDiscreteTraitTaxonSummaryRowImpl
    )
    assert (
        PosteriorMissingNucleotideDefinition is PosteriorMissingNucleotideDefinitionImpl
    )
    assert PosteriorMissingNucleotideReport is PosteriorMissingNucleotideReportImpl
    assert (
        PosteriorMissingNucleotideSequenceRecord
        is PosteriorMissingNucleotideSequenceRecordImpl
    )
    assert (
        PosteriorMissingNucleotideSiteSummaryRow
        is PosteriorMissingNucleotideSiteSummaryRowImpl
    )
    assert (
        PosteriorMissingNucleotideStateProbabilityRow
        is PosteriorMissingNucleotideStateProbabilityRowImpl
    )
    assert (
        build_posterior_missing_continuous_trait_definition
        is build_posterior_missing_continuous_trait_definition_impl
    )
    assert (
        build_posterior_missing_discrete_trait_definition
        is build_posterior_missing_discrete_trait_definition_impl
    )
    assert (
        build_posterior_missing_nucleotide_definition
        is build_posterior_missing_nucleotide_definition_impl
    )
    assert (
        summarize_brownian_continuous_trait_posterior_missing_values
        is summarize_brownian_continuous_trait_posterior_missing_values_impl
    )
    assert (
        summarize_continuous_trait_posterior_missing_values
        is summarize_continuous_trait_posterior_missing_values_impl
    )
    assert (
        summarize_discrete_trait_mk_posterior_missing_states
        is summarize_discrete_trait_mk_posterior_missing_states_impl
    )
    assert (
        summarize_fixed_topology_dna_posterior_missing_states
        is summarize_fixed_topology_dna_posterior_missing_states_impl
    )
    assert (
        summarize_fixed_topology_partitioned_dna_posterior_missing_states
        is summarize_fixed_topology_partitioned_dna_posterior_missing_states_impl
    )
    assert (
        summarize_joint_topology_dna_posterior_missing_states
        is summarize_joint_topology_dna_posterior_missing_states_impl
    )
    assert (
        summarize_nucleotide_posterior_missing_states
        is summarize_nucleotide_posterior_missing_states_impl
    )
    assert (
        summarize_ornstein_uhlenbeck_continuous_trait_posterior_missing_values
        is summarize_ornstein_uhlenbeck_continuous_trait_posterior_missing_values_impl
    )
