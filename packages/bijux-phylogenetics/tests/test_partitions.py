from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.core.partitions import (
    parse_locus_partitions,
    slice_partition_sequence,
    validate_locus_partitions,
    write_locus_partitions,
)
from bijux_phylogenetics.errors import InvalidPartitionError


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
