from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian.partition_model_priors import (
    build_partition_substitution_model_definition,
)
from bijux_phylogenetics.bayesian.posterior_ancestral_sequences import (
    PosteriorAncestralSequenceDefinition,
    build_posterior_ancestral_sequence_definition,
)
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    LocusSegment,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def _alignment_records() -> tuple[AlignmentRecord, ...]:
    return (
        AlignmentRecord(identifier="A", sequence="ACGT"),
        AlignmentRecord(identifier="B", sequence="ACGT"),
    )


def _locus_partition(name: str, start: int, end: int) -> LocusPartition:
    return LocusPartition(
        name=name,
        segments=(LocusSegment(start=start, end=end),),
        total_sites=(end - start) + 1,
        data_type="DNA",
    )


def test_build_posterior_ancestral_sequence_definition_records_partition_surface() -> (
    None
):
    definition = build_posterior_ancestral_sequence_definition(
        records=_alignment_records(),
        posterior_probability_threshold=0.75,
        minimum_clade_posterior_support=0.6,
        low_confidence_state_symbol="?",
        locus_partitions=(
            _locus_partition("gene_alpha", 1, 2),
            _locus_partition("gene_beta", 3, 4),
        ),
        partition_models=(
            build_partition_substitution_model_definition(
                partition_name="gene_alpha",
                model_name="K80",
            ),
            build_partition_substitution_model_definition(
                partition_name="gene_beta",
                model_name="K80",
            ),
        ),
    )

    assert isinstance(definition, PosteriorAncestralSequenceDefinition)
    assert definition.posterior_probability_threshold == 0.75
    assert definition.minimum_clade_posterior_support == 0.6
    assert definition.low_confidence_state_symbol == "?"
    assert tuple(partition.name for partition in definition.locus_partitions or ()) == (
        "gene_alpha",
        "gene_beta",
    )


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        (
            {"records": ()},
            "at least one alignment record",
        ),
        (
            {
                "records": (
                    AlignmentRecord(identifier="A", sequence="ACGT"),
                    AlignmentRecord(identifier="B", sequence="ACG"),
                )
            },
            "equal sequence lengths",
        ),
        (
            {
                "records": _alignment_records(),
                "posterior_probability_threshold": 1.2,
            },
            "at most one",
        ),
        (
            {
                "records": _alignment_records(),
                "locus_partitions": (_locus_partition("gene_alpha", 1, 2),),
            },
            "requires locus partitions and partition models together",
        ),
        (
            {
                "records": _alignment_records(),
                "locus_partitions": (
                    _locus_partition("gene_alpha", 1, 2),
                    _locus_partition("gene_beta", 3, 3),
                ),
                "partition_models": (
                    build_partition_substitution_model_definition(
                        partition_name="gene_alpha",
                        model_name="K80",
                    ),
                    build_partition_substitution_model_definition(
                        partition_name="gene_beta",
                        model_name="K80",
                    ),
                ),
            },
            "cover every alignment site",
        ),
    ],
)
def test_build_posterior_ancestral_sequence_definition_rejects_invalid_inputs(
    kwargs: dict[str, object],
    match: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=match):
        build_posterior_ancestral_sequence_definition(**kwargs)
