from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.errors import InvalidAlignmentError
from bijux_phylogenetics.io.fasta import (
    inspect_coding_alignment,
    load_fasta_alignment,
    load_permissive_fasta_records,
    prepare_coding_sequences_for_alignment,
    summarise_fasta,
    translate_coding_alignment,
)
from bijux_phylogenetics.shared_dna_alignment_fixtures import (
    get_shared_dna_alignment_fixture,
    list_shared_dna_alignment_fixtures,
)


def test_shared_dna_alignment_fixture_catalog_covers_required_goal_cases() -> None:
    fixtures = list_shared_dna_alignment_fixtures()
    feature_tags = {tag for fixture in fixtures for tag in fixture.feature_tags}

    assert {
        "clean-aligned-dna",
        "gaps",
        "ambiguous-iupac",
        "lowercase-input",
        "identical-sequences",
        "high-divergence-sequences",
        "missing-data",
        "unequal-length-invalid-input",
        "valid-reading-frame",
        "frame-error",
        "internal-stop-codon",
        "terminal-stop-codon",
    } <= feature_tags


def test_shared_dna_alignment_fixture_lookup_preserves_durable_ids() -> None:
    fixture = get_shared_dna_alignment_fixture("coding_terminal_stop")

    assert (
        fixture.relative_path
        == "alignments/example_alignment_coding_terminal_stop.fasta"
    )
    assert "terminal-stop-codon" in fixture.feature_tags
    assert fixture.path.is_file()


@pytest.mark.parametrize(
    ("fixture_id", "expected_sequence_count", "expected_alignment_length"),
    [
        ("clean_aligned_dna", 4, 8),
        ("dna_with_gaps", 4, 6),
        ("dna_with_ambiguity", 3, 6),
        ("lowercase_aligned_dna", 3, 6),
        ("identical_sequences", 4, 8),
        ("high_divergence_sequences", 3, 4),
        ("dna_with_missing_data", 3, 6),
        ("coding_valid_reading_frame", 3, 9),
        ("coding_internal_stop", 1, 9),
        ("coding_terminal_stop", 1, 9),
    ],
)
def test_shared_dna_alignment_fixture_catalog_loads_valid_cases(
    fixture_id: str,
    expected_sequence_count: int,
    expected_alignment_length: int,
) -> None:
    fixture = get_shared_dna_alignment_fixture(fixture_id)

    records = load_fasta_alignment(fixture.path)
    summary = summarise_fasta(fixture.path)

    assert len(records) == expected_sequence_count
    assert summary.sequence_count == expected_sequence_count
    assert summary.alignment_length == expected_alignment_length
    assert summary.inferred_alphabet == "dna"


def test_shared_dna_alignment_fixture_catalog_preserves_lowercase_records() -> None:
    fixture = get_shared_dna_alignment_fixture("lowercase_aligned_dna")

    records = load_fasta_alignment(fixture.path)

    assert records[0].sequence == "acgtaa"
    assert records[1].sequence == "acgtta"


def test_shared_dna_alignment_fixture_catalog_reports_unequal_length_invalid_input() -> (
    None
):
    fixture = get_shared_dna_alignment_fixture("unequal_length_invalid_input")

    records = load_permissive_fasta_records(fixture.path)

    assert [len(record.sequence) for record in records] == [6, 5]
    with pytest.raises(InvalidAlignmentError):
        load_fasta_alignment(fixture.path)


def test_shared_dna_alignment_fixture_catalog_translates_valid_reading_frame() -> None:
    fixture = get_shared_dna_alignment_fixture("coding_valid_reading_frame")

    translated, report = translate_coding_alignment(fixture.path)

    assert [record.sequence for record in translated] == ["MEL", "MKM", "MTW"]
    assert report.translated_sequence_count == 3
    assert report.stop_codon_count == 0
    prepared, preparation_report = prepare_coding_sequences_for_alignment(fixture.path)
    assert [record.sequence for record in prepared] == [
        "ATGGAACTG",
        "ATGAAAATG",
        "ATGACCTGG",
    ]
    assert preparation_report.warnings == []


def test_shared_dna_alignment_fixture_catalog_rejects_frame_error_translation() -> None:
    fixture = get_shared_dna_alignment_fixture("coding_frame_error")

    with pytest.raises(InvalidAlignmentError):
        translate_coding_alignment(fixture.path)
    with pytest.raises(InvalidAlignmentError):
        prepare_coding_sequences_for_alignment(fixture.path)

    diagnostics = inspect_coding_alignment(fixture.path)
    assert len(diagnostics.frameshift_like_sequences) == 1
    assert diagnostics.coding_behaviors[0].identifier == "frame_error"


def test_shared_dna_alignment_fixture_catalog_marks_internal_stop_behavior() -> None:
    fixture = get_shared_dna_alignment_fixture("coding_internal_stop")

    translated, report = translate_coding_alignment(fixture.path)

    assert [record.sequence for record in translated] == ["M*W"]
    assert report.stop_codon_count == 1
    diagnostics = inspect_coding_alignment(fixture.path)
    assert diagnostics.coding_behaviors[0].premature_stop_count == 1
    with pytest.raises(InvalidAlignmentError):
        prepare_coding_sequences_for_alignment(fixture.path)


def test_shared_dna_alignment_fixture_catalog_retains_terminal_stop_behavior() -> None:
    fixture = get_shared_dna_alignment_fixture("coding_terminal_stop")

    translated, report = translate_coding_alignment(fixture.path)

    assert [record.sequence for record in translated] == ["ME*"]
    assert report.stop_codon_count == 1
    diagnostics = inspect_coding_alignment(fixture.path)
    assert diagnostics.coding_behaviors[0].terminal_stop_count == 1
    prepared, preparation_report = prepare_coding_sequences_for_alignment(fixture.path)
    assert [record.sequence for record in prepared] == ["ATGGAATAA"]
    assert preparation_report.terminal_stop_sequence_count == 1
    assert (
        "terminal stop codons were retained in accepted coding sequences"
        in preparation_report.warnings
    )
