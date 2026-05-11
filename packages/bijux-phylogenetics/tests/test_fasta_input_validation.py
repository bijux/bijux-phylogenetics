from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.core.alignment import AlignmentRecord
from bijux_phylogenetics.errors import InvalidAlignmentError
from bijux_phylogenetics.io.fasta import (
    load_fasta_records,
    load_permissive_fasta_records,
    repair_fasta_input,
    validate_fasta_input,
    write_fasta_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(path: str) -> Path:
    return FIXTURES / path


def test_permissive_loader_preserves_duplicate_and_empty_raw_records() -> None:
    records = load_permissive_fasta_records(
        fixture("alignments/example_sequences_invalid_input.fasta")
    )

    assert [(record.identifier, record.sequence) for record in records] == [
        ("Alpha sample", "ACTGACTG"),
        ("Alpha sample", "ACTG1CTG"),
        ("odd id", ""),
        ("rare taxon", "ACTGACTGACTGACTGACTGACTG"),
    ]


def test_strict_loader_rejects_duplicate_identifiers_in_raw_input() -> None:
    with pytest.raises(InvalidAlignmentError, match="duplicate sequence ids"):
        load_fasta_records(fixture("alignments/example_sequences_invalid_input.fasta"))


def test_validate_fasta_input_reports_real_input_problems() -> None:
    report = validate_fasta_input(
        fixture("alignments/example_sequences_invalid_input.fasta"),
        sequence_type="dna",
    )

    assert report.summary.sequence_count == 4
    assert report.summary.empty_sequence_count == 1
    assert report.summary.min_sequence_length == 0
    assert report.summary.max_sequence_length == 24
    assert [row.identifier for row in report.duplicate_identifiers] == ["Alpha sample"]
    assert [
        (row.identifier, row.position, row.character)
        for row in report.illegal_characters
    ] == [("Alpha sample", 5, "1")]
    assert [(row.identifier, row.record_index) for row in report.empty_sequences] == [
        ("odd id", 3)
    ]
    assert [row.identifier for row in report.length_outliers] == ["odd id", "rare taxon"]
    assert report.warnings == [
        "input contains duplicate sequence identifiers",
        "input contains unsupported sequence characters",
        "input contains empty sequences",
        "input contains sequence length outliers",
    ]


def test_repair_fasta_input_normalizes_ids_and_removes_invalid_records() -> None:
    repaired_records, report = repair_fasta_input(
        fixture("alignments/example_sequences_invalid_input.fasta"),
        sequence_type="dna",
        normalize_identifiers=True,
        remove_invalid_records=True,
    )

    assert [(record.identifier, record.sequence) for record in repaired_records] == [
        ("Alpha_sample", "ACTGACTG"),
        ("rare_taxon", "ACTGACTGACTGACTGACTGACTG"),
    ]
    assert report.before.sequence_count == 4
    assert report.after.sequence_count == 2
    assert [row.repaired_identifier for row in report.normalized_identifiers] == [
        "Alpha_sample",
        "rare_taxon",
    ]
    assert [(row.identifier, row.reason) for row in report.removed_records] == [
        ("Alpha sample", "illegal-characters"),
        ("odd id", "empty-sequence"),
    ]


def test_repair_fasta_input_resolves_duplicate_collisions_when_ids_are_normalized(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "collisions.fasta"
    write_fasta_alignment(
        input_path,
        [
            AlignmentRecord(identifier="Alpha sample", sequence="ACTG"),
            AlignmentRecord(identifier="Alpha/sample", sequence="ACTA"),
            AlignmentRecord(identifier="Beta", sequence="ACTC"),
        ],
    )

    repaired_records, report = repair_fasta_input(
        input_path,
        sequence_type="dna",
        normalize_identifiers=True,
        remove_invalid_records=False,
    )

    assert [record.identifier for record in repaired_records] == [
        "Alpha_sample",
        "Alpha_sample_2",
        "Beta",
    ]
    assert [row.note for row in report.normalized_identifiers] == [
        "normalized identifier",
        "normalized identifier and resolved duplicate collision",
    ]
