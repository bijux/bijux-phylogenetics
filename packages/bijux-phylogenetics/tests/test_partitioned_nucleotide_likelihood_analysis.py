from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    LocusSegment,
    slice_partition_sequence,
)
from bijux_phylogenetics.phylo.likelihood import (
    PartitionedNucleotideModelAssignment,
    compile_partitioned_nucleotide_likelihood_analysis,
    evaluate_compiled_finite_state_likelihood_analysis,
    evaluate_compiled_finite_state_partition_log_likelihood_rows,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    validate_compiled_finite_state_partition_likelihood_reconstruction,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_partitioned_nucleotide_analysis_reuses_linked_hky_context() -> None:
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )
    partition_map = _two_locus_partition_map()
    assignments = (
        PartitionedNucleotideModelAssignment(
            partition_name="locus-alpha",
            model_name="hky85",
            context_linkage_key="shared-hky",
            kappa=4.0,
            base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
        ),
        PartitionedNucleotideModelAssignment(
            partition_name="locus-beta",
            model_name="hky85",
            context_linkage_key="shared-hky",
            kappa=4.0,
            base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
        ),
    )
    analysis = compile_partitioned_nucleotide_likelihood_analysis(
        tree,
        records,
        partition_map=partition_map,
        partition_model_assignments=assignments,
    )

    rows, row_total = evaluate_compiled_finite_state_partition_log_likelihood_rows(
        analysis
    )
    expected_total = math.fsum(
        evaluate_hky85_tree_likelihood(
            tree,
            _slice_partition_records(records, partition),
            kappa=4.0,
            base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
        ).log_likelihood
        for partition in partition_map
    )

    assert len(analysis.contexts) == 1
    assert [row.context_linkage_key for row in rows] == ["shared-hky", "shared-hky"]
    assert [row.context_index for row in rows] == [0, 0]
    validate_compiled_finite_state_partition_likelihood_reconstruction(
        rows,
        expected_total_log_likelihood=row_total,
        expected_partition_count=analysis.partition_count,
        owner_name="partitioned nucleotide linked hky test",
    )
    assert math.isclose(row_total, expected_total, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        row_total,
        evaluate_compiled_finite_state_likelihood_analysis(analysis),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_partitioned_nucleotide_analysis_keeps_unlinked_models_distinct() -> None:
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )
    partition_map = _two_locus_partition_map()
    analysis = compile_partitioned_nucleotide_likelihood_analysis(
        tree,
        records,
        partition_map=partition_map,
        partition_model_assignments=(
            PartitionedNucleotideModelAssignment(
                partition_name="locus-alpha",
                model_name="jc69",
            ),
            PartitionedNucleotideModelAssignment(
                partition_name="locus-beta",
                model_name="hky85",
                kappa=4.0,
                base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
            ),
        ),
    )

    rows, row_total = evaluate_compiled_finite_state_partition_log_likelihood_rows(
        analysis
    )
    expected_total = (
        evaluate_jc69_tree_likelihood(
            tree,
            _slice_partition_records(records, partition_map[0]),
        ).log_likelihood
        + evaluate_hky85_tree_likelihood(
            tree,
            _slice_partition_records(records, partition_map[1]),
            kappa=4.0,
            base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
        ).log_likelihood
    )

    assert len(analysis.contexts) == 2
    assert [row.partition_name for row in rows] == ["locus-alpha", "locus-beta"]
    assert [row.context_linkage_key for row in rows] == ["locus-alpha", "locus-beta"]
    assert rows[0].model_name != rows[1].model_name
    assert math.isclose(row_total, expected_total, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        row_total,
        evaluate_compiled_finite_state_likelihood_analysis(analysis),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_partitioned_nucleotide_analysis_requires_full_assignment_coverage() -> None:
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )

    with pytest.raises(
        ValueError,
        match="requires one model assignment for every declared partition",
    ):
        compile_partitioned_nucleotide_likelihood_analysis(
            tree,
            records,
            partition_map=_two_locus_partition_map(),
            partition_model_assignments=(
                PartitionedNucleotideModelAssignment(
                    partition_name="locus-alpha",
                    model_name="jc69",
                ),
            ),
        )


def _two_locus_partition_map() -> tuple[LocusPartition, ...]:
    return (
        LocusPartition(
            name="locus-alpha",
            segments=(LocusSegment(start=1, end=2),),
            total_sites=2,
            data_type="DNA",
        ),
        LocusPartition(
            name="locus-beta",
            segments=(LocusSegment(start=3, end=4),),
            total_sites=2,
            data_type="DNA",
        ),
    )


def _slice_partition_records(
    records: list[AlignmentRecord],
    partition: LocusPartition,
) -> list[AlignmentRecord]:
    return [
        AlignmentRecord(
            identifier=record.identifier,
            sequence=slice_partition_sequence(record.sequence, partition),
        )
        for record in records
    ]
