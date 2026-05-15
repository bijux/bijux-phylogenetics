from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.distance import (
    compute_pairwise_genetic_distance_matrix,
    compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment,
)
from bijux_phylogenetics.errors import InvalidAlignmentError
from bijux_phylogenetics.io.fasta import (
    compute_alignment_base_frequency_report,
    compute_alignment_base_frequency_report_from_dna_bin_alignment,
    compute_alignment_segregating_site_report,
    compute_alignment_segregating_site_report_from_dna_bin_alignment,
    load_dna_bin_alignment,
    write_dna_bin_alignment_fasta,
)


def fixture(name: str) -> Path:
    return Path(__file__).resolve().parent / "fixtures" / "alignments" / name


def test_load_dna_bin_alignment_preserves_taxon_order_length_and_case_normalization() -> (
    None
):
    alignment = load_dna_bin_alignment(fixture("example_alignment_lowercase.fasta"))

    assert alignment.source_alphabet == "dna"
    assert alignment.sequence_count == 3
    assert alignment.alignment_length == 6
    assert [record.identifier for record in alignment.records] == [
        "lower_a",
        "lower_b",
        "lower_c",
    ]
    assert [record.sequence for record in alignment.records] == [
        "acgtaa",
        "acgtta",
        "atgtnn",
    ]


def test_load_dna_bin_alignment_preserves_gap_missing_and_ambiguity_states() -> None:
    alignment = load_dna_bin_alignment(fixture("example_alignment_ambiguity.fasta"))

    assert [record.sequence for record in alignment.records] == [
        "acgtn?",
        "acgtr?",
        "acgt-?",
    ]
    assert [row.state for row in alignment.rows if row.identifier == "B"] == [
        "a",
        "c",
        "g",
        "t",
        "r",
        "?",
    ]


def test_load_dna_bin_alignment_rejects_invalid_symbols_explicitly() -> None:
    with pytest.raises(
        InvalidAlignmentError,
        match="dnabin-compatible nucleotide loading: A:5=Z",
    ):
        load_dna_bin_alignment(fixture("example_alignment_invalid_dna.fasta"))


def test_write_dna_bin_alignment_fasta_roundtrips_without_state_loss(
    tmp_path: Path,
) -> None:
    alignment = load_dna_bin_alignment(fixture("example_alignment_ambiguity.fasta"))
    output_path = tmp_path / "roundtrip.fasta"

    write_dna_bin_alignment_fasta(output_path, alignment)
    reloaded = load_dna_bin_alignment(output_path)

    assert reloaded.records == alignment.records
    assert reloaded.rows == alignment.rows


def test_dna_bin_alignment_supports_base_frequency_report_without_reloading() -> None:
    alignment = load_dna_bin_alignment(fixture("example_alignment_lowercase.fasta"))

    direct_report = compute_alignment_base_frequency_report_from_dna_bin_alignment(
        alignment
    )
    path_report = compute_alignment_base_frequency_report(
        fixture("example_alignment_lowercase.fasta")
    )

    assert direct_report.path == alignment.path
    assert direct_report.sequence_count == alignment.sequence_count
    assert direct_report.alignment_length == alignment.alignment_length
    assert direct_report.alignment_rows == path_report.alignment_rows
    assert direct_report.per_sequence_rows == path_report.per_sequence_rows
    assert direct_report.composition_outliers == path_report.composition_outliers


def test_dna_bin_alignment_supports_segregating_site_report_without_reloading() -> None:
    alignment = load_dna_bin_alignment(fixture("example_alignment_ambiguity.fasta"))

    direct_report = compute_alignment_segregating_site_report_from_dna_bin_alignment(
        alignment
    )
    path_report = compute_alignment_segregating_site_report(
        fixture("example_alignment_ambiguity.fasta")
    )

    assert direct_report.path == alignment.path
    assert direct_report.sequence_count == alignment.sequence_count
    assert direct_report.alignment_length == alignment.alignment_length
    assert direct_report.segregating_site_positions == path_report.segregating_site_positions
    assert direct_report.rows == path_report.rows
    assert direct_report.warnings == path_report.warnings


@pytest.mark.parametrize(
    "model",
    ["raw", "jc69", "k80", "f81", "tn93"],
)
def test_dna_bin_alignment_supports_nucleotide_distance_models_without_reloading(
    model: str,
) -> None:
    alignment = load_dna_bin_alignment(fixture("example_alignment_distance_gaps.fasta"))

    direct_report = compute_pairwise_genetic_distance_matrix_from_dna_bin_alignment(
        alignment,
        model=model,
    )
    path_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_gaps.fasta"),
        model=model,
    )

    assert direct_report.path == alignment.path
    assert direct_report.identifiers == [record.identifier for record in alignment.records]
    assert direct_report.alignment_length == alignment.alignment_length
    assert direct_report.model == path_report.model
    assert direct_report.model_parameters == path_report.model_parameters
    assert direct_report.warnings == path_report.warnings
    assert direct_report.pairs == path_report.pairs


@pytest.mark.parametrize(
    "model",
    ["raw", "jc69", "k80", "f81", "tn93"],
)
def test_all_dna_distance_models_reject_invalid_symbols_consistently(
    model: str,
) -> None:
    with pytest.raises(
        InvalidAlignmentError,
        match="dnabin-compatible nucleotide loading: A:5=Z",
    ):
        compute_pairwise_genetic_distance_matrix(
            fixture("example_alignment_invalid_dna.fasta"),
            model=model,
        )
