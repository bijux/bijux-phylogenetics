from __future__ import annotations

from bijux_phylogenetics.bayesian import (
    PosteriorAncestralSequenceDefinition,
    PosteriorAncestralSequenceRecord,
    PosteriorAncestralSequenceReport,
    PosteriorAncestralSiteSummaryRow,
    PosteriorAncestralStateProbabilityRow,
    build_posterior_ancestral_sequence_definition,
    summarize_nucleotide_posterior_ancestral_sequences,
    write_posterior_ancestral_sequence_fasta,
    write_posterior_ancestral_state_probability_table,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    PosteriorAncestralSequenceDefinition as PosteriorAncestralSequenceDefinitionImpl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    PosteriorAncestralSequenceRecord as PosteriorAncestralSequenceRecordImpl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    PosteriorAncestralSequenceReport as PosteriorAncestralSequenceReportImpl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    PosteriorAncestralSiteSummaryRow as PosteriorAncestralSiteSummaryRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    PosteriorAncestralStateProbabilityRow as PosteriorAncestralStateProbabilityRowImpl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    build_posterior_ancestral_sequence_definition as build_posterior_ancestral_sequence_definition_impl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    summarize_nucleotide_posterior_ancestral_sequences as summarize_nucleotide_posterior_ancestral_sequences_impl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    write_posterior_ancestral_sequence_fasta as write_posterior_ancestral_sequence_fasta_impl,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    write_posterior_ancestral_state_probability_table as write_posterior_ancestral_state_probability_table_impl,
)


def test_bayesian_exports_posterior_ancestral_sequence_surface() -> None:
    assert (
        PosteriorAncestralSequenceDefinition is PosteriorAncestralSequenceDefinitionImpl
    )
    assert PosteriorAncestralSequenceRecord is PosteriorAncestralSequenceRecordImpl
    assert PosteriorAncestralSequenceReport is PosteriorAncestralSequenceReportImpl
    assert PosteriorAncestralSiteSummaryRow is PosteriorAncestralSiteSummaryRowImpl
    assert (
        PosteriorAncestralStateProbabilityRow
        is PosteriorAncestralStateProbabilityRowImpl
    )
    assert (
        build_posterior_ancestral_sequence_definition
        is build_posterior_ancestral_sequence_definition_impl
    )
    assert (
        summarize_nucleotide_posterior_ancestral_sequences
        is summarize_nucleotide_posterior_ancestral_sequences_impl
    )
    assert (
        write_posterior_ancestral_sequence_fasta
        is write_posterior_ancestral_sequence_fasta_impl
    )
    assert (
        write_posterior_ancestral_state_probability_table
        is write_posterior_ancestral_state_probability_table_impl
    )
