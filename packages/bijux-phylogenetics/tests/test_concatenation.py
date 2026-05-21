from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from bijux_phylogenetics.io.fasta import write_fasta_alignment
import bijux_phylogenetics.phylo.alignment.concatenation as concatenation_api
from bijux_phylogenetics.phylo.alignment.concatenation import (
    concatenate_locus_alignments,
)
from bijux_phylogenetics.phylo.alignment.partitions import write_locus_partitions
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def _matrix_tsv(report: Any) -> str:
    def render(value: object) -> str:
        if isinstance(value, bool):
            return str(value).lower()
        if isinstance(value, float):
            return format(value, ".12g")
        return str(value)

    locus_names = [partition.name for partition in report.partitions]
    lines = [
        "\t".join(
            [
                "taxon",
                *locus_names,
                "covered_locus_count",
                "total_locus_count",
                "locus_coverage_fraction",
                "observed_site_count",
                "total_site_count",
                "site_coverage_fraction",
                "low_coverage",
            ]
        )
    ]
    lines.extend(
        "\t".join(
            render(value)
            for value in [
                row.taxon,
                *[row.occupancies[name] for name in locus_names],
                row.covered_locus_count,
                row.total_locus_count,
                row.locus_coverage_fraction,
                row.observed_site_count,
                row.total_site_count,
                row.site_coverage_fraction,
                row.low_coverage,
            ]
        )
        for row in report.taxa
    )
    return "\n".join(lines) + "\n"


def test_concatenate_locus_alignments_matches_expected_supermatrix_outputs(
    tmp_path: Path,
) -> None:
    input_dir = fixture("concatenation/mixed-locus-supermatrix/inputs")
    alignment_paths = [
        input_dir / "alpha-dna.fasta",
        input_dir / "beta-protein.fasta",
        input_dir / "gamma-dna.fasta",
    ]

    records, partitions, report = concatenate_locus_alignments(
        alignment_paths,
        data_types=("DNA", "PROTEIN", "DNA"),
        concatenated_alignment_path=tmp_path / "mixed-locus-supermatrix.aln.fasta",
        concatenated_partition_path=tmp_path / "mixed-locus-supermatrix.partitions.txt",
    )

    alignment_path = tmp_path / "mixed-locus-supermatrix.aln.fasta"
    partition_path = tmp_path / "mixed-locus-supermatrix.partitions.txt"
    write_fasta_alignment(alignment_path, records)
    write_locus_partitions(partition_path, partitions)

    assert alignment_path.read_text(encoding="utf-8") == fixture(
        "expected/concatenation/mixed-locus-supermatrix.aln.fasta"
    ).read_text(encoding="utf-8")
    assert partition_path.read_text(encoding="utf-8") == fixture(
        "expected/concatenation/mixed-locus-supermatrix.partitions.txt"
    ).read_text(encoding="utf-8")
    assert _matrix_tsv(report.occupancy_report) == fixture(
        "expected/concatenation/mixed-locus-supermatrix.matrix.tsv"
    ).read_text(encoding="utf-8")


def test_concatenate_locus_alignments_reports_missing_taxa_and_mixed_data_types() -> (
    None
):
    input_dir = fixture("concatenation/mixed-locus-supermatrix/inputs")
    _, partitions, report = concatenate_locus_alignments(
        [
            input_dir / "alpha-dna.fasta",
            input_dir / "beta-protein.fasta",
            input_dir / "gamma-dna.fasta",
        ],
        data_types=("DNA", "PROTEIN", "DNA"),
    )

    assert report.taxa == ["TaxonA", "TaxonB", "TaxonC", "TaxonD", "TaxonE"]
    assert [partition.data_type for partition in partitions] == [
        "DNA",
        "PROTEIN",
        "DNA",
    ]
    assert report.loci[0].missing_taxa == ["TaxonD", "TaxonE"]
    assert report.loci[1].missing_taxa == ["TaxonB", "TaxonC"]
    assert report.loci[2].missing_taxa == ["TaxonC", "TaxonD", "TaxonE"]
    assert report.occupancy_report.taxa[1].occupancies == {
        "alpha-dna": 1.0,
        "beta-protein": 0.0,
        "gamma-dna": 1.0,
    }


def test_concatenate_locus_alignments_rejects_duplicate_taxa_within_one_locus(
    tmp_path: Path,
) -> None:
    locus_path = tmp_path / "alpha-dna.fasta"
    locus_path.write_text(
        ">TaxonA\nATGC\n>TaxonA\nATGT\n",
        encoding="utf-8",
    )

    with pytest.raises(
        InvalidAlignmentError,
        match="duplicate sequence ids: TaxonA",
    ):
        concatenate_locus_alignments([locus_path])


def test_package_root_exports_concatenation_workflow() -> None:
    assert (
        concatenation_api.concatenate_locus_alignments is concatenate_locus_alignments
    )
