from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.core.locus_occupancy import (
    build_locus_occupancy_report,
    parse_locus_partitions,
)
from bijux_phylogenetics.errors import InvalidPartitionError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def test_parse_locus_partitions_accepts_common_partition_syntax(tmp_path: Path) -> None:
    partition_path = tmp_path / "partitions.txt"
    partition_path.write_text(
        "charset gene_alpha = 1-4;\n"
        "DNA,gene_beta = 5-9\n"
        "gene_gamma = 10-12\n",
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
        InvalidPartitionError, match="partition name 'gene_alpha' appears more than once"
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
    assert taxon_rows["TaxonA"].locus_coverage_fraction == 1.0
    assert taxon_rows["TaxonB"].locus_coverage_fraction == 2 / 3
    assert taxon_rows["TaxonC"].locus_coverage_fraction == 1 / 3
    assert set(report.low_coverage_taxa) == {"TaxonC", "TaxonD", "TaxonE"}
    assert locus_rows["gene_alpha"].taxon_coverage_fraction == 0.6
    assert locus_rows["gene_beta"].taxon_coverage_fraction == 0.6
    assert locus_rows["gene_gamma"].taxon_coverage_fraction == 0.4
    assert report.low_coverage_loci == ["gene_gamma"]
    assert taxon_rows["TaxonA"].occupancies == {
        "gene_alpha": 1.0,
        "gene_beta": 1.0,
        "gene_gamma": 1.0,
    }
