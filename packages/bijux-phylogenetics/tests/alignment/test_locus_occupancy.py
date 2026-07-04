from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta import load_fasta_alignment, write_fasta_alignment
from bijux_phylogenetics.phylo.alignment.occupancy import (
    build_locus_occupancy_report,
    build_locus_occupancy_report_from_records,
    filter_locus_occupancy,
    parse_locus_partitions,
    write_locus_partitions,
)
from bijux_phylogenetics.runtime.errors import InvalidPartitionError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def test_parse_locus_partitions_accepts_common_partition_syntax(tmp_path: Path) -> None:
    partition_path = tmp_path / "partitions.txt"
    partition_path.write_text(
        "charset gene_alpha = 1-4;\nDNA,gene_beta = 5-9\ngene_gamma = 10-12\n",
        encoding="utf-8",
    )

    partitions = parse_locus_partitions(partition_path)

    assert [partition.name for partition in partitions] == [
        "gene_alpha",
        "gene_beta",
        "gene_gamma",
    ]
    assert [partition.total_sites for partition in partitions] == [4, 5, 3]
    assert partitions[1].data_type == "DNA"


def test_parse_locus_partitions_rejects_overlapping_ranges_in_report_build(
    tmp_path: Path,
) -> None:
    with pytest.raises(InvalidPartitionError):
        build_locus_occupancy_report(
            fixture("alignments/example_multilocus_alignment.fasta"),
            _write_overlapping_partition_file(tmp_path),
        )


def _write_overlapping_partition_file(tmp_path: Path) -> Path:
    path = tmp_path / "overlapping_partitions.txt"
    path.write_text("gene_alpha = 1-4\ngene_beta = 4-9\n", encoding="utf-8")
    return path


def test_parse_locus_partitions_rejects_duplicate_locus_names(
    tmp_path: Path,
) -> None:
    partition_path = tmp_path / "duplicate_names.txt"
    partition_path.write_text(
        "gene_alpha = 1-4\ngene_alpha = 5-9\n",
        encoding="utf-8",
    )

    with pytest.raises(
        InvalidPartitionError,
        match="partition name 'gene_alpha' appears more than once",
    ):
        build_locus_occupancy_report(
            fixture("alignments/example_multilocus_alignment.fasta"),
            partition_path,
        )


def test_build_locus_occupancy_report_tracks_taxon_and_locus_coverage() -> None:
    report = build_locus_occupancy_report(
        fixture("alignments/example_multilocus_alignment.fasta"),
        fixture("alignments/example_multilocus_partitions.txt"),
        taxon_coverage_threshold=0.5,
        locus_coverage_threshold=0.5,
    )

    taxon_rows = {row.taxon: row for row in report.taxa}
    locus_rows = {row.locus_name: row for row in report.loci}

    assert report.taxon_count == 5
    assert report.locus_count == 3
    assert report.assigned_site_count == 12
    assert report.unassigned_site_count == 0
    assert report.minimum_locus_occupancy == 0.0
    assert taxon_rows["TaxonA"].locus_coverage_fraction == 1.0
    assert taxon_rows["TaxonB"].site_coverage_fraction == 7 / 12
    assert taxon_rows["TaxonB"].locus_coverage_fraction == 2 / 3
    assert taxon_rows["TaxonC"].locus_coverage_fraction == 1 / 3
    assert set(report.low_coverage_taxa) == {"TaxonC", "TaxonD", "TaxonE"}
    assert locus_rows["gene_alpha"].taxon_coverage_fraction == 0.6
    assert locus_rows["gene_beta"].taxon_coverage_fraction == 0.6
    assert locus_rows["gene_gamma"].taxon_coverage_fraction == 0.4
    assert locus_rows["gene_gamma"].site_coverage_fraction == 0.4
    assert report.low_coverage_loci == ["gene_gamma"]
    assert taxon_rows["TaxonA"].occupancies == {
        "gene_alpha": 1.0,
        "gene_beta": 1.0,
        "gene_gamma": 1.0,
    }


def test_build_locus_occupancy_report_from_records_reuses_loaded_inputs() -> None:
    alignment_path = fixture("alignments/example_multilocus_alignment.fasta")
    partition_path = fixture("alignments/example_multilocus_partitions.txt")
    records = load_fasta_alignment(alignment_path)
    partitions = parse_locus_partitions(partition_path)

    report = build_locus_occupancy_report_from_records(
        records=records,
        partitions=partitions,
        alignment_path=alignment_path,
        partition_path=partition_path,
    )

    assert report.alignment_path == alignment_path
    assert report.partition_path == partition_path
    assert report.taxon_count == 5
    assert report.locus_count == 3


def test_build_locus_occupancy_report_supports_minimum_cell_occupancy() -> None:
    report = build_locus_occupancy_report(
        fixture("alignments/example_multilocus_partial_occupancy.fasta"),
        fixture("alignments/example_multilocus_partial_occupancy_partitions.txt"),
        taxon_coverage_threshold=0.5,
        locus_coverage_threshold=0.75,
        minimum_locus_occupancy=0.75,
    )

    taxon_rows = {row.taxon: row for row in report.taxa}
    locus_rows = {row.locus_name: row for row in report.loci}

    assert report.minimum_locus_occupancy == 0.75
    assert taxon_rows["TaxonB"].occupancies["gene_gamma"] == 0.25
    assert taxon_rows["TaxonB"].covered_locus_count == 2
    assert taxon_rows["TaxonB"].site_coverage_fraction == 0.75
    assert report.low_coverage_taxa == ["TaxonC", "TaxonD", "TaxonE"]
    assert locus_rows["gene_gamma"].covered_taxon_count == 2
    assert locus_rows["gene_gamma"].taxon_coverage_fraction == 0.4
    assert locus_rows["gene_gamma"].site_coverage_fraction == 0.45
    assert report.low_coverage_loci == ["gene_alpha", "gene_gamma"]


def test_filter_locus_occupancy_can_remove_low_coverage_loci_only() -> None:
    records, partitions, filter_report = filter_locus_occupancy(
        fixture("alignments/example_multilocus_alignment.fasta"),
        fixture("alignments/example_multilocus_partitions.txt"),
        locus_coverage_threshold=0.5,
    )

    assert [record.identifier for record in records] == [
        "TaxonA",
        "TaxonB",
        "TaxonC",
        "TaxonD",
        "TaxonE",
    ]
    assert [partition.name for partition in partitions] == ["gene_alpha", "gene_beta"]
    assert len(records[0].sequence) == 9
    assert filter_report.retained_loci == ["gene_alpha", "gene_beta"]
    assert filter_report.removed_loci == ["gene_gamma"]
    assert filter_report.final_report.locus_count == 2


def test_filter_locus_occupancy_iterates_until_taxa_and_loci_stabilize() -> None:
    records, partitions, filter_report = filter_locus_occupancy(
        fixture("alignments/example_multilocus_alignment.fasta"),
        fixture("alignments/example_multilocus_partitions.txt"),
        taxon_coverage_threshold=0.6,
        locus_coverage_threshold=0.6,
    )

    assert [record.identifier for record in records] == ["TaxonA", "TaxonB"]
    assert [record.sequence for record in records] == ["AAAAGGG", "AAAAGGG"]
    assert [partition.name for partition in partitions] == ["gene_alpha", "gene_gamma"]
    assert [partition.total_sites for partition in partitions] == [4, 3]
    assert filter_report.removed_taxa == ["TaxonC", "TaxonD", "TaxonE"]
    assert filter_report.removed_loci == ["gene_beta"]
    assert filter_report.filtered_alignment_length == 7
    assert filter_report.iterations == 2
    assert filter_report.final_report.low_coverage_taxa == []
    assert filter_report.final_report.low_coverage_loci == []


def test_filter_locus_occupancy_records_iteration_history() -> None:
    records, partitions, filter_report = filter_locus_occupancy(
        fixture("alignments/example_multilocus_partial_occupancy.fasta"),
        fixture("alignments/example_multilocus_partial_occupancy_partitions.txt"),
        taxon_coverage_threshold=0.5,
        locus_coverage_threshold=0.75,
        minimum_locus_occupancy=0.75,
    )

    assert [record.identifier for record in records] == ["TaxonA", "TaxonB"]
    assert [record.sequence for record in records] == ["AAAACCCC", "AAAACCCC"]
    assert [partition.name for partition in partitions] == ["gene_alpha", "gene_beta"]
    assert filter_report.minimum_locus_occupancy == 0.75
    assert filter_report.removed_taxa == ["TaxonC", "TaxonD", "TaxonE"]
    assert filter_report.removed_loci == ["gene_gamma"]
    assert filter_report.iterations == 2
    assert [step.iteration for step in filter_report.filter_iterations] == [1, 2]
    assert filter_report.filter_iterations[0].removed_taxa == [
        "TaxonC",
        "TaxonD",
        "TaxonE",
    ]
    assert filter_report.filter_iterations[0].removed_loci == ["gene_gamma"]
    assert filter_report.filter_iterations[1].removed_taxa == []
    assert filter_report.filter_iterations[1].removed_loci == []
    assert filter_report.final_report.low_coverage_taxa == []
    assert filter_report.final_report.low_coverage_loci == []


def test_locus_occupancy_rejects_invalid_minimum_cell_occupancy() -> None:
    with pytest.raises(
        ValueError,
        match="minimum locus occupancy must be between 0 and 1 inclusive",
    ):
        build_locus_occupancy_report(
            fixture("alignments/example_multilocus_alignment.fasta"),
            fixture("alignments/example_multilocus_partitions.txt"),
            minimum_locus_occupancy=1.1,
        )


def test_write_locus_partitions_persists_remapped_coordinates(tmp_path: Path) -> None:
    records, partitions, _ = filter_locus_occupancy(
        fixture("alignments/example_multilocus_alignment.fasta"),
        fixture("alignments/example_multilocus_partitions.txt"),
        taxon_coverage_threshold=0.6,
        locus_coverage_threshold=0.6,
    )
    alignment_path = tmp_path / "filtered_alignment.fasta"
    partition_path = tmp_path / "filtered_partitions.txt"

    write_fasta_alignment(alignment_path, records)
    write_locus_partitions(partition_path, partitions)

    assert alignment_path.read_text(encoding="utf-8") == (
        ">TaxonA\nAAAAGGG\n>TaxonB\nAAAAGGG\n"
    )
    assert partition_path.read_text(encoding="utf-8") == (
        "DNA,gene_alpha = 1-4\nDNA,gene_gamma = 5-7\n"
    )
