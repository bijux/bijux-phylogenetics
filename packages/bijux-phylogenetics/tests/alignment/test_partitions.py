from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.alignment.partitions import (
    build_partition_summary_report,
    parse_locus_partitions,
    slice_partition_sequence,
    validate_locus_partitions,
    write_locus_partitions,
    write_partition_summary_table,
)
from bijux_phylogenetics.runtime.errors import InvalidPartitionError


def test_parse_locus_partitions_accepts_whitespace_and_stride_syntax(
    tmp_path: Path,
) -> None:
    partition_path = tmp_path / "partitions.nex"
    partition_path.write_text(
        "#nexus\n"
        "begin sets;\n"
        "charset gene_alpha = 1-6\\3 2-6\\3;\n"
        "AA, gene_beta = 7-10, 12;\n"
        "end;\n",
        encoding="utf-8",
    )

    partitions = parse_locus_partitions(partition_path)

    assert [partition.name for partition in partitions] == ["gene_alpha", "gene_beta"]
    assert [partition.total_sites for partition in partitions] == [4, 5]
    assert partitions[0].segments[0].step == 3
    assert partitions[1].data_type == "PROTEIN"


def test_validate_locus_partitions_rejects_stride_overlap(tmp_path: Path) -> None:
    partition_path = tmp_path / "overlap.partitions"
    partition_path.write_text(
        "DNA,gene_alpha = 1-9\\3\nDNA,gene_beta = 4-6\n",
        encoding="utf-8",
    )

    partitions = parse_locus_partitions(partition_path)

    with pytest.raises(
        InvalidPartitionError,
        match="partition 'gene_beta' overlaps another locus at site 4",
    ):
        validate_locus_partitions(partitions, alignment_length=12)


def test_slice_partition_sequence_applies_segment_stride(tmp_path: Path) -> None:
    partition_path = tmp_path / "stride.partitions"
    partition_path.write_text("DNA,gene_alpha = 1-6\\3,2-6\\3\n", encoding="utf-8")
    partitions = parse_locus_partitions(partition_path)

    assert slice_partition_sequence("ABCDEFGHIJKL", partitions[0]) == "ADBE"


def test_write_locus_partitions_preserves_stride_syntax(tmp_path: Path) -> None:
    partition_path = tmp_path / "input.partitions"
    partition_path.write_text(
        "DNA,gene_alpha = 1-9\\3\nPROTEIN,gene_beta = 10-12\n",
        encoding="utf-8",
    )
    partitions = parse_locus_partitions(partition_path)
    output_path = tmp_path / "written.partitions"

    write_locus_partitions(output_path, partitions)

    assert output_path.read_text(encoding="utf-8") == (
        "DNA,gene_alpha = 1-9\\3\nPROTEIN,gene_beta = 10-12\n"
    )


def test_build_partition_summary_report_tracks_three_gene_partition_set(
    tmp_path: Path,
) -> None:
    partition_path = tmp_path / "summary.partitions"
    partition_path.write_text(
        "DNA,gene_alpha = 1-4\nDNA,gene_beta = 5-9\nPROTEIN,gene_gamma = 10-12\n",
        encoding="utf-8",
    )
    partitions = parse_locus_partitions(partition_path)

    report = build_partition_summary_report(partitions, alignment_length=14)

    assert report.partition_count == 3
    assert report.assigned_site_count == 12
    assert report.unassigned_site_count == 2
    assert report.mixed_data_types is True
    assert report.declared_data_types == ["DNA", "PROTEIN"]
    assert report.rows[1].coordinate_text == "5-9"
    assert report.warnings == ["2 alignment sites are not assigned to any partition"]


def test_write_partition_summary_table_persists_tsv_rows(tmp_path: Path) -> None:
    partition_path = tmp_path / "summary.partitions"
    partition_path.write_text(
        "DNA,gene_alpha = 1-4\nDNA,gene_beta = 5-9\nDNA,gene_gamma = 10-12\n",
        encoding="utf-8",
    )
    partitions = parse_locus_partitions(partition_path)
    report = build_partition_summary_report(partitions, alignment_length=12)
    output_path = tmp_path / "partition-summary.tsv"

    write_partition_summary_table(output_path, report)

    assert output_path.read_text(encoding="utf-8") == (
        "locus_name\tdata_type\tsegment_count\ttotal_sites\tstart_site\tend_site\tcoordinate_text\n"
        "gene_alpha\tDNA\t1\t4\t1\t4\t1-4\n"
        "gene_beta\tDNA\t1\t5\t5\t9\t5-9\n"
        "gene_gamma\tDNA\t1\t3\t10\t12\t10-12\n"
    )
