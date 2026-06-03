from __future__ import annotations

import pytest

from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_dna_alignment_fixture,
    list_shared_dna_alignment_fixtures,
)
from bijux_phylogenetics.io.fasta import (
    inspect_coding_alignment,
    load_dna_bin_alignment,
    load_fasta_alignment,
    load_permissive_fasta_records,
    prepare_coding_sequences_for_alignment,
    summarise_fasta,
    translate_coding_alignment,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError


def test_shared_dna_alignment_fixture_catalog_covers_required_goal_cases() -> None:
    fixtures = list_shared_dna_alignment_fixtures()
    feature_tags = {tag for fixture in fixtures for tag in fixture.feature_tags}

    assert {
        "clean-aligned-dna",
        "gaps",
        "ambiguous-iupac",
        "lowercase-input",
        "invariant-alignment",
        "one-variable-site",
        "identical-sequences",
        "high-divergence-sequences",
        "missing-data",
        "all-gap-missing",
        "invalid-symbol",
        "unequal-length-invalid-input",
        "valid-reading-frame",
        "frame-error",
        "internal-stop-codon",
        "ambiguous-codon",
        "terminal-stop-codon",
        "alternate-genetic-code",
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
    (
        "fixture_id",
        "expected_sequence_count",
        "expected_alignment_length",
        "expected_alphabet",
    ),
    [
        ("clean_aligned_dna", 4, 8, "dna"),
        ("dna_with_gaps", 4, 6, "dna"),
        ("dna_with_ambiguity", 3, 6, "dna"),
        ("lowercase_aligned_dna", 3, 6, "dna"),
        ("invariant_aligned_dna", 3, 6, "dna"),
        ("one_variable_site_alignment", 3, 6, "dna"),
        ("identical_sequences", 4, 8, "dna"),
        ("high_divergence_sequences", 3, 4, "dna"),
        ("dna_with_missing_data", 3, 6, "dna"),
        ("all_gap_missing_alignment", 3, 6, "unknown"),
        ("coding_valid_reading_frame", 3, 9, "dna"),
        ("coding_ambiguous_codon", 1, 9, "dna"),
        ("coding_internal_stop", 1, 9, "dna"),
        ("coding_terminal_stop", 1, 9, "dna"),
        ("coding_mitochondrial_triplet", 1, 9, "dna"),
    ],
)
def test_shared_dna_alignment_fixture_catalog_loads_valid_cases(
    fixture_id: str,
    expected_sequence_count: int,
    expected_alignment_length: int,
    expected_alphabet: str,
) -> None:
    fixture = get_shared_dna_alignment_fixture(fixture_id)

    records = load_fasta_alignment(fixture.path)
    summary = summarise_fasta(fixture.path)

    assert len(records) == expected_sequence_count
    assert summary.sequence_count == expected_sequence_count
    assert summary.alignment_length == expected_alignment_length
    assert summary.inferred_alphabet == expected_alphabet


def test_shared_dna_alignment_fixture_catalog_preserves_lowercase_records() -> None:
    fixture = get_shared_dna_alignment_fixture("lowercase_aligned_dna")

    records = load_fasta_alignment(fixture.path)
    matrix = load_dna_bin_alignment(fixture.path)

    assert records[0].sequence == "acgtaa"
    assert records[1].sequence == "acgtta"
    assert [record.sequence for record in matrix.records] == [
        "acgtaa",
        "acgtta",
        "atgtnn",
    ]


def test_shared_dna_alignment_fixture_catalog_marks_invalid_symbol_rejection() -> None:
    fixture = get_shared_dna_alignment_fixture("invalid_symbol_alignment")

    assert fixture.load_expectation == "invalid-symbol"
    with pytest.raises(InvalidAlignmentError, match="A:5=Z"):
        load_dna_bin_alignment(fixture.path)


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


def test_shared_dna_alignment_fixture_catalog_truncates_frame_error_translation() -> (
    None
):
    fixture = get_shared_dna_alignment_fixture("coding_frame_error")

    translated, report = translate_coding_alignment(fixture.path)

    assert [record.sequence for record in translated] == ["ME"]
    assert report.translated_alignment_length == 2
    assert report.dropped_trailing_nucleotide_count == 2
    assert report.trailing_partial_codon_sequence_count == 1
    assert report.warnings == [
        "sequence length not a multiple of 3: 2 nucleotides dropped"
    ]
    with pytest.raises(InvalidAlignmentError):
        prepare_coding_sequences_for_alignment(fixture.path)

    diagnostics = inspect_coding_alignment(fixture.path)
    assert len(diagnostics.frameshift_like_sequences) == 1
    assert diagnostics.coding_behaviors[0].identifier == "frame_error"


def test_shared_dna_alignment_fixture_catalog_marks_ambiguous_translation_behavior() -> (
    None
):
    fixture = get_shared_dna_alignment_fixture("coding_ambiguous_codon")

    translated, report = translate_coding_alignment(fixture.path)

    assert [record.sequence for record in translated] == ["MXG"]
    assert report.invalid_codon_count == 1
    assert report.stop_codon_count == 0
    assert (
        report.codon_observations[1].translation_status == "ambiguous-or-invalid-codon"
    )
    with pytest.raises(InvalidAlignmentError):
        prepare_coding_sequences_for_alignment(fixture.path)


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


def test_shared_dna_alignment_fixture_catalog_honors_alternate_genetic_code() -> None:
    fixture = get_shared_dna_alignment_fixture("coding_mitochondrial_triplet")

    standard, standard_report = translate_coding_alignment(
        fixture.path, genetic_code="1"
    )
    mitochondrial, mitochondrial_report = translate_coding_alignment(
        fixture.path,
        genetic_code="2",
    )

    assert [record.sequence for record in standard] == ["M*G"]
    assert [record.sequence for record in mitochondrial] == ["MWG"]
    assert standard_report.stop_codon_count == 1
    assert mitochondrial_report.stop_codon_count == 0
