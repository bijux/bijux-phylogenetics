from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta import (
    detect_fasta_sequence_type,
    load_fasta_records,
    load_permissive_fasta_records,
    repair_fasta_input,
    validate_fasta_input,
    write_fasta_alignment,
)
from bijux_phylogenetics.phylo.alignment import AlignmentRecord
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

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


def test_detect_fasta_sequence_type_classifies_supported_and_blocked_inputs(
    tmp_path: Path,
) -> None:
    dna_path = tmp_path / "dna.fasta"
    rna_path = tmp_path / "rna.fasta"
    protein_path = tmp_path / "protein.fasta"
    mixed_path = tmp_path / "mixed.fasta"

    write_fasta_alignment(
        dna_path,
        [
            AlignmentRecord(identifier="A", sequence="ACTGACTG"),
            AlignmentRecord(identifier="B", sequence="ACTTACTA"),
        ],
    )
    write_fasta_alignment(
        rna_path,
        [
            AlignmentRecord(identifier="A", sequence="ACUGACUG"),
            AlignmentRecord(identifier="B", sequence="ACUGACUA"),
        ],
    )
    write_fasta_alignment(
        protein_path,
        [
            AlignmentRecord(identifier="P1", sequence="MKTWFLIM"),
            AlignmentRecord(identifier="P2", sequence="MRTWYLVM"),
        ],
    )
    write_fasta_alignment(
        mixed_path,
        [
            AlignmentRecord(identifier="dna_like", sequence="ACTGACTG"),
            AlignmentRecord(identifier="rna_like", sequence="ACUGACUG"),
        ],
    )

    dna_report = detect_fasta_sequence_type(dna_path)
    assert dna_report.detected_type == "dna"
    assert dna_report.selected_type == "dna"
    assert dna_report.compatible_types == ["dna", "protein"]
    assert dna_report.confidence == "medium"

    rna_report = detect_fasta_sequence_type(rna_path)
    assert rna_report.detected_type == "rna"
    assert rna_report.selected_type == "rna"
    assert rna_report.compatible_types == ["rna"]
    assert rna_report.confidence == "high"

    protein_report = detect_fasta_sequence_type(protein_path)
    assert protein_report.detected_type == "protein"
    assert protein_report.selected_type == "protein"
    assert protein_report.compatible_types == ["protein"]
    assert protein_report.confidence == "high"

    mixed_report = detect_fasta_sequence_type(mixed_path)
    assert mixed_report.detected_type == "mixed"
    assert mixed_report.selected_type is None
    assert mixed_report.compatible_types == []
    assert mixed_report.confidence == "blocked"

    invalid_report = detect_fasta_sequence_type(
        fixture("alignments/example_sequences_invalid_input.fasta")
    )
    assert invalid_report.detected_type == "invalid"
    assert invalid_report.selected_type is None
    assert invalid_report.invalid_record_count == 1


def test_detect_fasta_sequence_type_handles_ambiguity_codes(tmp_path: Path) -> None:
    dna_path = tmp_path / "dna-ambiguity.fasta"
    rna_path = tmp_path / "rna-ambiguity.fasta"
    protein_path = tmp_path / "protein-ambiguity.fasta"

    write_fasta_alignment(
        dna_path,
        [
            AlignmentRecord(identifier="A", sequence="ACTGNRYM"),
            AlignmentRecord(identifier="B", sequence="ACTGBDHV"),
        ],
    )
    write_fasta_alignment(
        rna_path,
        [
            AlignmentRecord(identifier="A", sequence="ACUGNRYM"),
            AlignmentRecord(identifier="B", sequence="ACUGBDHV"),
        ],
    )
    write_fasta_alignment(
        protein_path,
        [
            AlignmentRecord(identifier="P1", sequence="MKTXBZX"),
            AlignmentRecord(identifier="P2", sequence="MRQXBZX"),
        ],
    )

    dna_report = detect_fasta_sequence_type(dna_path)
    assert dna_report.detected_type == "dna"
    assert dna_report.selected_type == "dna"

    rna_report = detect_fasta_sequence_type(rna_path)
    assert rna_report.detected_type == "rna"
    assert rna_report.selected_type == "rna"

    protein_report = detect_fasta_sequence_type(protein_path)
    assert protein_report.detected_type == "protein"
    assert protein_report.selected_type == "protein"


def test_validate_fasta_input_reports_real_input_problems() -> None:
    report = validate_fasta_input(
        fixture("alignments/example_sequences_invalid_input.fasta"),
        sequence_type="dna",
    )

    assert report.summary.sequence_count == 4
    assert report.sequence_type_report.detected_type == "invalid"
    assert report.sequence_type_report.selected_type is None
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
    assert [row.identifier for row in report.length_outliers] == [
        "odd id",
        "rare taxon",
    ]
    assert report.warnings == [
        "input contains duplicate sequence identifiers",
        "input contains unsupported sequence characters",
        "input contains empty sequences",
        "input contains sequence length outliers",
    ]


def test_validate_fasta_input_surfaces_automatic_sequence_type_warning() -> None:
    report = validate_fasta_input(fixture("alignments/example_sequences_raw.fasta"))

    assert report.summary.inferred_alphabet == "dna"
    assert report.sequence_type_report.detected_type == "dna"
    assert report.sequence_type_report.confidence == "medium"
    assert (
        "automatic sequence type defaults to dna from nucleotide-like characters "
        "that remain protein-compatible by alphabet alone"
    ) in report.warnings


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


def test_repair_fasta_input_requires_explicit_type_for_mixed_records(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "mixed.fasta"
    write_fasta_alignment(
        input_path,
        [
            AlignmentRecord(identifier="dna_like", sequence="ACTGACTG"),
            AlignmentRecord(identifier="rna_like", sequence="ACUGACUG"),
        ],
    )

    with pytest.raises(InvalidAlignmentError, match="explicit sequence_type"):
        repair_fasta_input(
            input_path,
            normalize_identifiers=False,
            remove_invalid_records=True,
        )


def test_repair_fasta_input_can_remove_records_incompatible_with_declared_type(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "mixed.fasta"
    write_fasta_alignment(
        input_path,
        [
            AlignmentRecord(identifier="dna_like", sequence="ACTGACTG"),
            AlignmentRecord(identifier="rna_like", sequence="ACUGACUG"),
            AlignmentRecord(identifier="protein_like", sequence="MKTWFLIM"),
        ],
    )

    repaired_records, report = repair_fasta_input(
        input_path,
        sequence_type="dna",
        normalize_identifiers=False,
        remove_invalid_records=True,
    )

    assert [record.identifier for record in repaired_records] == ["dna_like"]
    assert [(row.identifier, row.reason) for row in report.removed_records] == [
        ("rna_like", "illegal-characters+sequence-type-mismatch"),
        ("protein_like", "illegal-characters+sequence-type-mismatch"),
    ]
    assert (
        "repair removed records incompatible with the declared sequence type"
        in report.warnings
    )
